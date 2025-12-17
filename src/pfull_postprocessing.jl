module Postprocessing

import NCDatasets
# import ClimaDiagnostics.Writers as Writers
# import ..Writers: NetCDFWriter
import ..Writers
import ..Writers: NetCDFWriter, PfullCoordsStyle
import ClimaCore
import ClimaCore: Remapping, Geometry

"""
    write_cc_grid_to_regular_grid(
        writer::NetCDFWriter{CS},
    ) where {CS <: PfullCoordsStyle}

Given the NetCDF `writer` whose coordinate style is `PfullCoordsStyle`,
convert the all NetCDF files to a regular grid of longitude, latitude, and
pressure levels.
"""
function write_cc_grid_to_regular_grid(writer::NetCDFWriter)
    # TODO: I don't know why I can't put this in the type declaration
    writer.coordinates_style isa PfullCoordsStyle ||
        error("NetCDFWriter must have pressure coordinates style")

    _pfull_dir = joinpath(writer.output_dir, "_pfull_coords")
    # TODO: I am not sure how this works with restarts
    pfull_dir = joinpath(writer.output_dir, "pfull_coords")
    # If it already exists, then the pressure coordinates are converted already
    # TODO: Deal with restarts...
    # TODO: This can be relaxed a little bit by checking each file and seeing if
    # it exists or not.
    isdir(pfull_dir) && return nothing
    mkdir(pfull_dir)

    remapper = LevelRemapper(writer)
    if isdir(_pfull_dir)
        pfull_nc_filepaths = filter!(
            filepath -> isfile(filepath) && endswith(filepath, ".nc"),
            readdir(_pfull_dir, join = true),
        )
        for pfull_nc_filepath in pfull_nc_filepaths
            write_cc_grid_to_regular_grid(
                writer,
                remapper,
                pfull_dir,
                pfull_nc_filepath,
            )
        end
    end
end

"""
    write_cc_grid_to_regular_grid(
        writer::NetCDFWriter,
        remapper,
        output_dir,
        filepath::String,
    )

Write a NetCDF file consisting of a ClimaCore grid to a regular grid of
longitude, latitude, and pressure levels.
"""
function write_cc_grid_to_regular_grid(
    writer::NetCDFWriter,
    remapper,
    output_dir,
    filepath::String,
)
    isdir(filepath) && return nothing
    endswith(filepath, ".nc") || return nothing
    ds = NCDatasets.NCDataset(filepath)
    if !haskey(ds, "horizontal_index")
        close(ds)
        return nothing
    end

    times = Array(ds["time"])
    pfull_levels = Array(ds["pressure_level"])
    target_lon, target_lat = writer.hpts

    # Find variable name by finding the variable with 3 dimensions (time x pressure level x horizontal_index)
    not_dim_names = setdiff!(keys(ds), NCDatasets.dimnames(ds))
    varidx = findfirst(
        name ->
            NCDatasets.dimnames(ds[name]) ==
            ("time", "pressure_level", "horizontal_index"),
        not_dim_names,
    )
    varname = not_dim_names[varidx]

    # Do interpolation here!
    interpolated_data = interpolate(remapper, ds, varname)

    # Make new dataset
    # Anything that is not varname, pressure_levels, horizontal_index, lon, and lat
    # add again (which means that the last quantity is time)
    filename = basename(filepath)
    nc_hcoords_filepath = joinpath(output_dir, filename)
    nc_hcoords = NCDatasets.NCDataset(nc_hcoords_filepath, "c")

    # Add dimensions (lon, lat, pressure_level) and its attributes
    # maybe use add_dimension function?
    attribs_to_nt = attrib -> NamedTuple(Symbol(k) => v for (k, v) in attrib)
    time_attribs = attribs_to_nt(ds["time"].attrib)

    # TODO: Import these...
    Writers.add_dimension!(nc_hcoords, "time", times; time_attribs...)
    time_bnds_attribs = attribs_to_nt(ds["time_bnds"].attrib)
    Writers.add_time_bounds_maybe!(
        nc_hcoords,
        eltype(times);
        time_bnds_attribs...,
    )
    nc_hcoords["time_bnds"][:, :] = Array(ds["time_bnds"][:, :])

    if "date" in keys(ds)
        dates = Array(ds["date"])

        NCDatasets.defDim(nc_hcoords, "date", size(dates)[end])

        date_dim = NCDatasets.defVar(
            nc_hcoords,
            "date",
            Float64,
            ("date",),
            attrib = ds["date"].attrib,
        )

        date_dim[:] = dates

        date_bnds_attribs = attribs_to_nt(ds["date_bnds"].attrib)
        Writers.add_date_bounds_maybe!(nc_hcoords; date_bnds_attribs...)
        nc_hcoords["date_bnds"][:, :] = Array(ds["date_bnds"][:, :])
    end

    Writers.add_dimension!(
        nc_hcoords,
        "lon",
        target_lon;
        units = "degrees_east",
        axis = "X",
        standard_name = "longitude",
        long_name = "Longitude",
    )
    Writers.add_dimension!(
        nc_hcoords,
        "lat",
        target_lat;
        units = "degrees_north",
        axis = "Y",
        standard_name = "latitude",
        long_name = "Latitude",
    )

    pfull_attribs = attribs_to_nt(ds["pressure_level"].attrib)
    Writers.add_dimension!(
        nc_hcoords,
        "pressure_level",
        pfull_levels;
        pfull_attribs...,
    )
    # Add variable of interest and its attributes
    v = NCDatasets.defVar(
        nc_hcoords,
        varname,
        eltype(interpolated_data),
        ("time", "pressure_level", "lon", "lat"),
        deflatelevel = writer.compression_level,
        attrib = ds[varname].attrib,
    )
    v[:, :, :, :] = interpolated_data

    close(ds)
    close(nc_hcoords)
    return nothing
end

# Do not need to store times because times is not necessarily the same between
# all of them
# Can store horizontal points because it is the same across all NetCDFWriter
struct LevelRemapper{R, FL <: ClimaCore.Fields.Field, HPTS}
    """ClimaCore remapper object defined on a level of the pressure field."""
    cc_remapper::R

    """A level of the pressure field."""
    field_level::FL

    """The horizontal points to interpolate to."""
    hpts::HPTS
end

function LevelRemapper(writer::NetCDFWriter)
    pfull_field = writer.coordinates_style.pressure_field
    pfull_field_level = ClimaCore.to_cpu(ClimaCore.level(pfull_field, 1))

    # Create Remapper object with the specified horizontal and vertical points
    space = axes(pfull_field_level)
    lons, lats = writer.hpts
    target_hcoords =
        [Geometry.LatLongPoint(lat, lon) for lon in lons, lat in lats]
    remapper = Remapping.Remapper(space, target_hcoords)
    return LevelRemapper(remapper, pfull_field_level, writer.hpts)
end

function interpolate(
    remapper::LevelRemapper,
    ds::NCDatasets.NCDataset,
    varname::String,
)
    data = Array(ds[varname])
    dimnames = NCDatasets.dimnames(ds[varname])
    return interpolate(remapper, data, dimnames)
end

function interpolate(remapper::LevelRemapper, data::Array, dimnames)
    ndims(data) == length(dimnames) || error(
        "Number of dimensions of data is not the same as the length of dimnames",
    )
    dim_lengths = Dict(zip(dimnames, size(data)))
    lons, lats = remapper.hpts
    lon_length = length(lons)
    lat_length = length(lats)
    # TODO: Names are hardcoded
    time_length = dim_lengths["time"]
    pfull_length = dim_lengths["pressure_level"]
    interpolated_data = zeros(time_length, pfull_length, lon_length, lat_length)

    # TODO: Ordering of the dimensions are hardcoded
    field_level = remapper.field_level
    for (idx, level) in pairs(eachslice(data, dims = (1, 2)))
        values = ClimaCore.DataLayouts.array2data(
            level,
            ClimaCore.Fields.field_values(field_level),
        )
        field = ClimaCore.Fields.Field(values, axes(field_level))
        dest = @view interpolated_data[idx, :, :]
        Remapping.interpolate!(dest, remapper.cc_remapper, field)
    end
    return interpolated_data
end

end
