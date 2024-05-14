"""
    add_dimension!(nc::NCDatasets.NCDataset,
                   name::String,
                   points;
                   kwargs...)

Add dimension identified by `name` in the given `nc` file and fill it with the given
`points`.
"""
function add_dimension!(
    nc::NCDatasets.NCDataset,
    name::String,
    points;
    kwargs...,
)
    FT = eltype(points)

    NCDatasets.defDim(nc, name, size(points)[end])

    dim = NCDatasets.defVar(nc, name, FT, (name,))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end

    dim[:] = points

    return nothing
end

"""
    dimension_exists(
        nc::NCDatasets.NCDataset,
        name::String,
        expected_size::Tuple,
        )

Return whether the given dimension exists in the given dataset, and if yes, it has the same
size as `expected_size`.
"""
function dimension_exists(
    nc::NCDatasets.NCDataset,
    name::String,
    expected_size::Tuple,
)
    if haskey(nc, name)
        if size(nc[name]) != expected_size
            file_path = NCDatasets.path(nc)
            error(
                "Incompatible $name dimension already exists in file $file_path",
            )
        else
            return true
        end
    else
        return false
    end
end

"""
    add_time_maybe!(nc::NCDatasets.NCDataset,
                    float_type::Type{FT};
                    kwargs...) where {FT}

Add the `time` dimension (with infinite size) to the given NetCDF file if not already there.
Optionally, add all the keyword arguments as attributes.

Also add a `date` dataset (as a string).
"""
function add_time_maybe!(
    nc::NCDatasets.NCDataset,
    float_type::Type{FT};
    kwargs...,
) where {FT}

    # If we already have time, do nothing
    haskey(nc, "time") && return nothing

    NCDatasets.defDim(nc, "time", Inf)
    dim = NCDatasets.defVar(nc, "time", FT, ("time",))
    NCDatasets.defVar(nc, "date", String, ("time",))
    for (k, v) in kwargs
        dim.attrib[String(k)] = v
    end
    return nothing
end

"""
    add_space_coordinates_maybe!(nc::NCDatasets.NCDataset,
                                 space::Spaces.AbstractSpace,
                                 num_points;
                                 names)

Add dimensions relevant to the `space` to the given `nc` NetCDF file. The range is
automatically determined and the number of points is set with `num_points`, which has to be
an iterable of size N, where N is the number of dimensions of the space. For instance, 3 for
a cubed sphere, 2 for a surface, 1 for a column.

The function returns an array with the names of the relevant dimensions. (We want arrays
because we want to preserve the order to match the one in num_points).

In some cases, the names are adjustable passing the keyword `names`.
"""
function add_space_coordinates_maybe! end

"""
    target_coordinates!(space::Spaces.AbstractSpace,
                        num_points)

Return the range of interpolation coordinates. The range is automatically determined and the
number of points is set with `num_points`, which has to be an iterable of size N, where N is
the number of dimensions of the space. For instance, 3 for a cubed sphere, 2 for a surface,
1 for a column.
"""
function target_coordinates(space, num_points) end

function target_coordinates(
    space::S,
    num_points;
    disable_vertical_interpolation,
) where {
    S <:
    Union{Spaces.CenterFiniteDifferenceSpace, Spaces.FaceFiniteDifferenceSpace},
}
    
    # HACK 1
    if disable_vertical_interpolation
        cspace = Spaces.space(space, Grids.CellCenter())
        return Array(parent(Fields.coordinate_field(cspace).z))[:, 1]
    end
    # 

    # Exponentially spaced with base e
    #
    # We mimic something that looks like pressure levels
    #
    # p ~ p₀ exp(-z/H)
    #
    # We assume H to be 7000, which is a good scale height for the Earth atmosphere
    H_EARTH = 7000

    num_points_z = num_points[]

    # HACK 2
    FT = Spaces.undertype(space)
    topology = Spaces.topology(space)
    vert_domain = topology.mesh.domain
    z_min, z_max = FT(vert_domain.coord_min.z), FT(vert_domain.coord_max.z)
    # We floor z_min to avoid having to deal with the singular value z = 0.
    z_min = max(z_min, 100)
    exp_z_min = exp(-z_min / H_EARTH)
    exp_z_max = exp(-z_max / H_EARTH)

    return collect(-H_EARTH * log.(range(exp_z_min, exp_z_max, num_points_z)))
end

# Column
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.FiniteDifferenceSpace,
    num_points_z;
    disable_vertical_interpolation,
    names = ("z",),
)
    name, _... = names
    z_dimension_exists = dimension_exists(nc, name, (num_points_z,))
    if !z_dimension_exists
        zpts = target_coordinates(
            space,
            num_points_z;
            disable_vertical_interpolation,
        )
        add_dimension!(nc, name, zpts, units = "m", axis = "Z")
    end
    return [name]
end

add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.AbstractSpectralElementSpace,
    num_points;
) = add_space_coordinates_maybe!(
    nc,
    space,
    num_points,
    Meshes.domain(Spaces.topology(space));
)


# For the horizontal space, we also have to look at the domain, so we define another set of
# functions that dispatches over the domain
target_coordinates(space::Spaces.AbstractSpectralElementSpace, num_points) =
    target_coordinates(space, num_points, Meshes.domain(Spaces.topology(space)))

# Box
function target_coordinates(
    space::Spaces.SpectralElementSpace2D,
    num_points,
    domain::Domains.RectangleDomain,
)
    num_points_x, num_points_y = num_points
    FT = Spaces.undertype(space)
    xmin = FT(domain.interval1.coord_min.x)
    xmax = FT(domain.interval1.coord_max.x)
    ymin = FT(domain.interval2.coord_min.y)
    ymax = FT(domain.interval2.coord_max.y)
    xpts = collect(range(xmin, xmax, num_points_x))
    ypts = collect(range(ymin, ymax, num_points_y))
    return (xpts, ypts)
end

# Plane
function target_coordinates(
    space::Spaces.SpectralElementSpace1D,
    num_points,
    domain::Domains.IntervalDomain,
)
    num_points_x, _... = num_points
    FT = Spaces.undertype(space)
    xmin = FT(domain.coord_min.x)
    xmax = FT(domain.coord_max.x)
    xpts = collect(range(xmin, xmax, num_points_x))
    return (xpts)
end

# Cubed sphere
function target_coordinates(
    space::Spaces.SpectralElementSpace2D,
    num_points,
    ::Domains.SphereDomain,
)
    num_points_long, num_points_lat = num_points
    FT = Spaces.undertype(space)
    longpts = collect(range(FT(-180), FT(180), num_points_long))
    latpts = collect(range(FT(-90), FT(90), num_points_lat))

    return (longpts, latpts)
end

# Box
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace2D,
    num_points,
    ::Domains.RectangleDomain;
    names = ("x", "y"),
)
    xname, yname = names
    num_points_x, num_points_y = num_points
    x_dimension_exists = dimension_exists(nc, xname, (num_points_x,))
    y_dimension_exists = dimension_exists(nc, yname, (num_points_y,))

    if !x_dimension_exists && !y_dimension_exists
        xpts, ypts = target_coordinates(space, num_points)
        add_dimension!(nc, xname, xpts; units = "m", axis = "X")
        add_dimension!(nc, yname, ypts; units = "m", axis = "Y")
    end
    return [xname, yname]
end

# Plane
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace1D,
    num_points,
    ::Domains.IntervalDomain;
    names = ("x",),
)
    xname, _... = names
    num_points_x, = num_points
    x_dimension_exists = dimension_exists(nc, xname, (num_points_x,))

    if !x_dimension_exists
        xpts = target_coordinates(space, num_points)
        add_dimension!(nc, xname, xpts; units = "m", axis = "X")
    end
    return [xname]
end

# Cubed sphere
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.SpectralElementSpace2D,
    num_points,
    ::Domains.SphereDomain;
    names = ("lon", "lat"),
)
    longname, latname = names
    num_points_long, num_points_lat = num_points

    long_dimension_exists = dimension_exists(nc, longname, (num_points_long,))
    lat_dimension_exists = dimension_exists(nc, latname, (num_points_lat,))

    if !long_dimension_exists && !lat_dimension_exists
        longpts, latpts = target_coordinates(space, num_points)
        add_dimension!(
            nc,
            longname,
            longpts;
            units = "degrees_east",
            axis = "X",
        )
        add_dimension!(nc, latname, latpts; units = "degrees_north", axis = "Y")
    end

    return [longname, latname]
end

# General hybrid space. This calls both the vertical and horizontal add_space_coordinates_maybe!
# and combines the resulting dictionaries
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.ExtrudedFiniteDifferenceSpace,
    num_points;
    disable_vertical_interpolation,
    interpolated_physical_z = nothing,
)

    hdims_names = vdims_names = []

    num_points_horiz..., num_points_vertic = num_points

    # Being an Extruded space, we can assume that we have an horizontal and a vertical space.
    # We can also assume that the vertical space has dimension 1
    horizontal_space = Spaces.horizontal_space(space)

    hdims_names =
        add_space_coordinates_maybe!(nc, horizontal_space, num_points_horiz)

    vertical_space = Spaces.FiniteDifferenceSpace(
        Spaces.vertical_topology(space),
        Spaces.staggering(space),
    )

    if Spaces.grid(space).hypsography isa Grids.Flat
        vdims_names = add_space_coordinates_maybe!(
            nc,
            vertical_space,
            num_points_vertic;
            disable_vertical_interpolation,
        )
    else
        vdims_names = add_space_coordinates_maybe!(
            nc,
            vertical_space,
            num_points_vertic,
            interpolated_physical_z;
            disable_vertical_interpolation,
            names = ("z_reference",),
            depending_on_dimensions = hdims_names,
        )
    end

    return (hdims_names..., vdims_names...)
end

# Ignore the interpolated_physical_z/disable_vertical_interpolation keywords in the general
# case (we only case about the specialized one for extruded spaces)
add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space,
    num_points;
    interpolated_physical_z = nothing,
    disable_vertical_interpolation = false,
) = add_space_coordinates_maybe!(nc::NCDatasets.NCDataset, space, num_points)

# Elevation with topography

# `depending_on_dimensions` identifies the dimensions upon which the current one depends on
# (excluding itself). In pretty much all cases, the dimensions depend only on themselves
# (e.g., `lat` is a variable only defined on the latitudes.), and `depending_on_dimensions`
# should be an empty tuple. The only case in which this is not what happens is with `z` with
# topography. With topography, the altitude will depend on the spatial coordinates. So,
# `depending_on_dimensions` might be `("lon", "lat)`, or similar.
function add_space_coordinates_maybe!(
    nc::NCDatasets.NCDataset,
    space::Spaces.FiniteDifferenceSpace,
    num_points,
    interpolated_physical_z;
    names = ("z_reference",),
    disable_vertical_interpolation,
    depending_on_dimensions,
)
    num_points_z = num_points
    name, _... = names

    # Add z_reference
    z_reference_dimension_dimension_exists =
        dimension_exists(nc, name, (num_points_z,))

    if !z_reference_dimension_dimension_exists
        reference_altitudes = target_coordinates(
            space,
            num_points_z;
            disable_vertical_interpolation,
        )
        add_dimension!(nc, name, reference_altitudes; units = "m", axis = "Z")
    end

    # We also have to add an extra variable with the physical altitudes
    physical_name = "z_physical"
    z_physical_dimension_dimension_exists =
        dimension_exists(nc, physical_name, size(interpolated_physical_z))

    if !z_physical_dimension_dimension_exists
        FT = eltype(interpolated_physical_z)
        dim = NCDatasets.defVar(
            nc,
            physical_name,
            FT,
            (depending_on_dimensions..., name),
        )
        dim.attrib["units"] = "m"
        if length(depending_on_dimensions) == 2
            dim[:, :, :] = interpolated_physical_z
        elseif length(depending_on_dimensions) == 1
            dim[:, :] = interpolated_physical_z
        else
            error("Error in calculating z_physical")
        end
    end
    # We do not output this name because it is not an axis

    return [name]
end

# General hybrid space. This calls both the vertical and horizontal add_space_coordinates_maybe!
# and combines the resulting dictionaries
function target_coordinates(
    space::Spaces.ExtrudedFiniteDifferenceSpace,
    num_points;
    disable_vertical_interpolation,
)

    hcoords = vcoords = ()

    num_points_horiz..., num_points_vertic = num_points

    hcoords =
        target_coordinates(Spaces.horizontal_space(space), num_points_horiz)

    vertical_space = Spaces.FiniteDifferenceSpace(
        Spaces.vertical_topology(space),
        Spaces.staggering(space),
    )
    vcoords = target_coordinates(
        vertical_space,
        num_points_vertic;
        disable_vertical_interpolation,
    )

    hcoords == vcoords == () && error("Found empty space")

    return hcoords, vcoords
end

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace2D,
    domain::Domains.SphereDomain,
    hpts,
)
    # Notice LatLong not LongLat!
    return [Geometry.LatLongPoint(hc2, hc1) for hc1 in hpts[1], hc2 in hpts[2]]
end

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace2D,
    domain::Domains.RectangleDomain,
    hpts,
)
    return [Geometry.XYPoint(hc1, hc2) for hc1 in hpts[1], hc2 in hpts[2]]
end

function hcoords_from_horizontal_space(
    space::Spaces.SpectralElementSpace1D,
    domain::Domains.IntervalDomain,
    hpts,
)
    return [Geometry.XPoint(hc1) for hc1 in hpts]
end

"""
    hcoords_from_horizontal_space(space, domain, hpts)

Prepare the matrix of horizontal coordinates with the correct type according to the given `space`
and `domain` (e.g., `ClimaCore.Geometry.LatLongPoint`s).
"""
function hcoords_from_horizontal_space(space, domain, hpts) end
