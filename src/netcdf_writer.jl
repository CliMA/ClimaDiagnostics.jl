import Dates

import ClimaCore: Domains, Geometry, Grids, Fields, Meshes, Spaces
import ClimaCore.Remapping: Remapper, interpolate, interpolate!

import NCDatasets

# Defines target_coordinates, add_space_coordinates_maybe!, add_time_maybe! for a bunch of
# Spaces
include("netcdf_writer_coordinates.jl")

"""
    NetCDFWriter

A struct to remap `ClimaCore` `Fields` to rectangular grids and save the output to NetCDF
files.
"""
struct NetCDFWriter{T, TS, DI} <: AbstractWriter
    """The base folder where to save the files."""
    output_dir::String

    # TODO: At the moment, each variable gets its remapper. This is a little bit of a waste
    # because we probably only need a handful of remappers since the same remapper can be
    # used for multiple fields as long as they are all defined on the same space. We need
    # just a few remappers because realistically we need to support fields defined on the
    # entire space and fields defined on 2D slices. However, handling this distinction at
    # construction time is quite difficult.
    """ClimaCore `Remapper`s that interpolate Fields to rectangular grids."""
    remappers::Dict{String, Remapper}

    """ Tuple/Array of integers that identifies how many points to use for interpolation
    along the various dimensions. It has to have the same size as the target interpolation
    space."""
    num_points::T

    """How much to compress the data in the final NetCDF file: 0 no compression, 9 max
    compression."""
    compression_level::Int

    """An array with size num_points with the physical altitude of any given target
    point."""
    interpolated_physical_z::TS

    """NetCDF files that are currently open. Only the root process uses this field."""
    open_files::Dict{String, NCDatasets.NCDataset}

    """ Do not interpolate on the z direction, instead evaluate on the levels. When
    disable_vertical_interpolation is true, the num_points on the vertical direction is
    ignored."""
    disable_vertical_interpolation::Bool

    """Areas of memory preallocated where the interpolation output is saved. Only the root
    process uses this."""
    preallocated_output_arrays::DI
end

"""
    close(writer::NetCDFWriter)

Close all the open files in `writer`.
"""
function Base.close(writer::NetCDFWriter)
    foreach(NCDatasets.close, values(writer.open_files))
    return nothing
end

"""
    NetCDFWriter(cspace, output_dir)

Save a `ScheduledDiagnostic` to a NetCDF file inside the `output_dir` of the simulation by
performing a pointwise (non-conservative) remapping first.

Keyword arguments
==================

- `space`: `Space` where the `Fields` are defined.
- `output_dir`: The base folder where the files should be saved.
- `num_points`: How many points to use along the different dimensions to interpolate the
                fields. This is a tuple of integers, typically having meaning Long-Lat-Z,
                or X-Y-Z (the details depend on the configuration being simulated).
- `disable_vertical_interpolation`: Do not interpolate on the z direction, instead evaluate
                                    at on levels. When disable_vertical_interpolation is true,
                                    the num_points on the vertical direction is ignored.
- `compression_level`: How much to compress the output NetCDF file (0 is no compression, 9
  is maximum compression).

"""
function NetCDFWriter(
    space,
    output_dir;
    num_points = (180, 90, 50),
    disable_vertical_interpolation = false,
    compression_level = 9,
)
    horizontal_space = Spaces.horizontal_space(space)
    is_horizontal_space = horizontal_space == space

    if disable_vertical_interpolation
        # It is a little tricky to override the number of vertical points because we don't
        # know if the vertical direction is the 2nd (as in a plane) or 3rd index (as in a
        # box or sphere). To set this value, we check if we are on a plane or not

        # TODO: Get the number of dimensions directly from the space
        num_horiz_dimensions =
            Spaces.horizontal_space(space) isa Spaces.SpectralElementSpace1D ?
            1 : 2

        num_vpts = Meshes.nelements(Grids.vertical_topology(space).mesh)

        @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
        num_points = Tuple([num_points[1:num_horiz_dimensions]..., num_vpts])
    end

    # Interpolate physical zs
    if is_horizontal_space
        hpts = target_coordinates(space, num_points)
        vpts = []
    else
        hpts, vpts = target_coordinates(
            space,
            num_points;
            disable_vertical_interpolation,
        )
    end

    hcoords = hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    zcoords = Geometry.ZPoint.(vpts)

    remapper = Remapper(space, hcoords, zcoords)

    interpolated_physical_z =
        interpolate(remapper, Fields.coordinate_field(space).z)

    preallocated_arrays =
        ClimaComms.iamroot(ClimaComms.context(space)) ? Dict{String, Array}() :
        Dict{String, Nothing}()

    return NetCDFWriter{
        typeof(num_points),
        typeof(interpolated_physical_z),
        typeof(preallocated_arrays),
    }(
        output_dir,
        Dict{String, Remapper}(),
        num_points,
        compression_level,
        interpolated_physical_z,
        Dict{String, NCDatasets.NCDataset}(),
        disable_vertical_interpolation,
        preallocated_arrays,
    )
end

function interpolate_field!(writer::NetCDFWriter, field, diagnostic, u, p, t)

    var = diagnostic.variable

    space = axes(field)

    horizontal_space = Spaces.horizontal_space(space)

    # We have to deal with to cases: when we have an horizontal slice (e.g., the
    # surface), and when we have a full space. We distinguish these cases by checking if
    # the given space has the horizontal_space attribute. If not, it is going to be a
    # SpectralElementSpace2D and we don't have to deal with the z coordinates.
    is_horizontal_space = horizontal_space == space

    # Prepare the remapper if we don't have one for the given variable. We need one remapper
    # per variable (not one per diagnostic since all the time reductions return the same
    # type of space).

    # TODO: Expand this once we support spatial reductions
    if !haskey(writer.remappers, var.short_name)

        # hpts, vpts are ranges of numbers
        # hcoords, zcoords are ranges of Geometry.Points

        zcoords = []

        if is_horizontal_space
            hpts = target_coordinates(space, writer.num_points)
            vpts = []
        else
            hpts, vpts = target_coordinates(
                space,
                writer.num_points;
                writer.disable_vertical_interpolation,
            )
        end

        hcoords = hcoords_from_horizontal_space(
            horizontal_space,
            Meshes.domain(Spaces.topology(horizontal_space)),
            hpts,
        )

        # When we disable vertical_interpolation, we override the vertical points with
        # the reference values for the vertical space.
        if writer.disable_vertical_interpolation && !is_horizontal_space
            # We need Array(parent()) because we want an array of values, not a DataLayout
            # of Points
            vpts = Array(
                parent(
                    space.grid.vertical_grid.center_local_geometry.coordinates,
                ),
            )
        end

        zcoords = [Geometry.ZPoint(p) for p in vpts]

        writer.remappers[var.short_name] = Remapper(space, hcoords, zcoords)
    end

    remapper = writer.remappers[var.short_name]

    # Now we can interpolate onto the target points
    # There's an MPI call in here (to aggregate the results)
    #
    # The first time we call this, we call interpolate and allocate a new array.
    # Future calls are in-place
    if haskey(writer.preallocated_output_arrays, var.short_name)
        interpolate!(
            writer.preallocated_output_arrays[var.short_name],
            remapper,
            field,
        )
    else
        writer.preallocated_output_arrays[var.short_name] =
            interpolate(remapper, field)
    end
    return nothing
end

"""
    write_field!(writer::NetCDFWriter, field::Fields.Field, diagnostic, u, p, t)

Save the resampled `field` produced by `diagnostic` as directed by the `writer`.

Only the root process does something here.

Note: It assumes that the field is already resampled.

The target file is determined by `output_short_name(diagnostic)`. If the target file already
exists, append to it. If not, create a new file. If the file does not contain dimensions,
they are added the first time something is written.

Attributes are appended to the dataset:
- `short_name`
- `long_name`
- `units`
- `comments`
- `start_date`
"""
function write_field!(writer::NetCDFWriter, field, diagnostic, u, p, t)
    # Only the root process has to write
    ClimaComms.iamroot(ClimaComms.context(field)) || return nothing

    var = diagnostic.variable
    interpolated_field = writer.preallocated_output_arrays[var.short_name]
    space = axes(field)
    FT = Spaces.undertype(space)

    output_path =
        joinpath(writer.output_dir, "$(output_short_name(diagnostic)).nc")

    if !haskey(writer.open_files, output_path)
        # Append or write a new file
        open_mode = isfile(output_path) ? "a" : "c"
        writer.open_files[output_path] =
            NCDatasets.Dataset(output_path, open_mode)
    end

    nc = writer.open_files[output_path]

    # Define time coordinate
    add_time_maybe!(nc, FT; units = "s", axis = "T")

    dim_names = add_space_coordinates_maybe!(
        nc,
        space,
        writer.num_points;
        writer.disable_vertical_interpolation,
        writer.interpolated_physical_z,
    )

    if haskey(nc, "$(var.short_name)")
        # We already have something in the file
        v = nc["$(var.short_name)"]
        temporal_size, spatial_size... = size(v)
        spatial_size == size(interpolated_field) ||
            error("incompatible dimensions for $(var.short_name)")
    else
        v = NCDatasets.defVar(
            nc,
            "$(var.short_name)",
            FT,
            ("time", dim_names...),
            deflatelevel = writer.compression_level,
        )
        v.attrib["short_name"] = var.short_name::String
        v.attrib["long_name"] = output_long_name(diagnostic)::String
        v.attrib["units"] = var.units::String
        v.attrib["comments"] = var.comments::String
        # FIXME: We are hardcoding p.start_date !
        v.attrib["start_date"] = string(p.start_date)::String
        temporal_size = 0
    end

    # We need to write to the next position after what we read from the data (or the first
    # position ever if we are writing the file for the first time)
    time_index = temporal_size + 1

    nc["time"][time_index] = t

    # FIXME: We are hardcoding p.start_date !
    # FIXME: We are rounding t
    nc["date"][time_index] = string(p.start_date + Dates.Second(round(t)))

    # TODO: It would be nice to find a cleaner way to do this
    if length(dim_names) == 3
        v[time_index, :, :, :] = interpolated_field
    elseif length(dim_names) == 2
        v[time_index, :, :] = interpolated_field
    elseif length(dim_names) == 1
        v[time_index, :] = interpolated_field
    end
    return nothing
end

function Base.show(io::IO, writer::NetCDFWriter)
    num_open_files = length(keys(writer.open_files))
    print(
        io,
        "NetCDFWriter, writing to $(writer.output_dir) ($num_open_files files open)",
    )
end