module Postprocessing

import NCDatasets
import ClimaDiagnostics.Writers as Writers
import ..Writers: NetCDFWriter
import ClimaCore
import ClimaCore: Remapping, Geometry
import ..ClimaDiagnostics: PfullCoordsDiagnosticsHandler

# TODO: Would be better if this was a string instead of nc?
# Make a new netcdf file with lon, lat, pfull (TODO: Check if this is the same as ERA5)
# function write_h_indices_to_h_coords(
#     filepath::String,
#     target_lon,
#     target_lat;
#     output_dir = nothing,
# )
#     isdir(filepath) && return nothing
#     endswith(filepath, ".nc") || return nothing
#     ds = NCDatasets.NCDataset(filepath)
#     if !haskey(ds, "horizontal_index")
#         close(ds)
#         return nothing
#     end

#     # Interpolate to regular grid using constant interpolation
#     times = Array(ds["time"])
#     lon_of_cols = Array(ds["lon"])
#     lat_of_cols = Array(ds["lat"])
#     pfull_levels = Array(ds["pressure_level"])

#     col_indices =
#         find_closest_col.(
#             target_lon,
#             target_lat',
#             Ref(lon_of_cols),
#             Ref(lat_of_cols),
#         )

#     # Find variable name by finding the variable with 3 dimensions (time x pressure level x horizontal_index)
#     not_dim_names = setdiff!(keys(ds), NCDatasets.dimnames(ds))
#     varidx = findfirst(
#         name ->
#             NCDatasets.dimnames(ds[name]) ==
#             ("time", "pressure_level", "horizontal_index"),
#         not_dim_names,
#     )
#     varname = not_dim_names[varidx]

#     data = Array(ds[varname])

#     # Constant interpolation along the horizontal direction
#     interpolated_data = data[:, :, col_indices]

#     # Make new dataset
#     # Anything that is not varname, pressure_levels, horizontal_index, lon, and lat
#     # add again (which means that the last quantity is time)
#     filename = basename(filepath)
#     nc_hcoords_filepath = joinpath(output_dir, filename)
#     nc_hcoords = NCDatasets.NCDataset(nc_hcoords_filepath, "c")

#     # Add dimensions (lon, lat, pressure_level) and its attributes
#     # maybe use add_dimension function?
#     attribs_to_nt = attrib -> NamedTuple(Symbol(k) => v for (k, v) in attrib)
#     time_attribs = attribs_to_nt(ds["time"].attrib)
#     add_dimension!(nc_hcoords, "time", times; time_attribs...)
#     time_bnds_attribs = attribs_to_nt(ds["time_bnds"].attrib)
#     add_time_bounds_maybe!(nc_hcoords, eltype(times); time_bnds_attribs...)
#     nc_hcoords["time_bnds"][:, :] = Array(ds["time_bnds"][:, :])

#     if "date" in keys(ds)
#         dates = Array(ds["date"])
#         # date_attribs = attribs_to_nt(ds["date"].attrib)
#         # add_dimension!(nc_hcoords, "date", dates; date_attribs...)

#         NCDatasets.defDim(nc_hcoords, "date", size(dates)[end])

#         date_dim = NCDatasets.defVar(
#             nc_hcoords,
#             "date",
#             Float64,
#             ("date",),
#             attrib = ds["date"].attrib,
#         )

#         date_dim[:] = dates

#         date_bnds_attribs = attribs_to_nt(ds["date_bnds"].attrib)
#         add_date_bounds_maybe!(nc_hcoords; date_bnds_attribs...)
#         nc_hcoords["date_bnds"][:, :] = Array(ds["date_bnds"][:, :])
#     end

#     add_dimension!(
#         nc_hcoords,
#         "lon",
#         target_lon;
#         units = "degrees_east",
#         axis = "X",
#         standard_name = "longitude",
#         long_name = "Longitude",
#     )
#     add_dimension!(
#         nc_hcoords,
#         "lat",
#         target_lat;
#         units = "degrees_north",
#         axis = "Y",
#         standard_name = "latitude",
#         long_name = "Latitude",
#     )

#     pfull_attribs = attribs_to_nt(ds["pressure_level"].attrib)
#     add_dimension!(nc_hcoords, "pressure_level", pfull_levels; pfull_attribs...)

#     # Add variable of interest and its attributes
#     v = NCDatasets.defVar(
#         nc_hcoords,
#         varname,
#         eltype(interpolated_data),
#         ("time", "pressure_level", "lon", "lat"),
#         # deflatelevel = writer.compression_level,
#         attrib = ds[varname].attrib,
#     )
#     v[:, :, :, :] = interpolated_data

#     close(ds)
#     close(nc_hcoords)
#     return nothing
# end

# # TODO: Ask Dennis about this, since I am not sure if there is a better way of doing it
# function find_closest_col(target_lon, target_lat, lon_of_cols, lat_of_cols)
#     # TODO: Maybe account for the fact that this is on a sphere
#     dist(lon, lat) = (lon - target_lon)^2 + (lat - target_lat)^2
#     min_idx = 1
#     min_dist = dist(first(lon_of_cols), first(lat_of_cols))
#     for (i, (lon, lat)) in enumerate(zip(lon_of_cols, lat_of_cols))
#         curr_dist = dist(lon, lat)
#         if curr_dist < min_dist
#             min_idx = i
#             min_dist = curr_dist
#         end
#     end
#     return min_idx
# end

###############################

function write_h_indices_to_regular_grid(
    pfull_diag_handler::PfullCoordsDiagnosticsHandler,
    writer::NetCDFWriter,
)
    remapper = create_pfull_coords_remapper(pfull_diag_handler, writer)
    _pfull_dir = joinpath(writer.output_dir, "_pfull_coords")
    # TODO: I am not sure how this works with restarts
    pfull_dir = joinpath(writer.output_dir, "pfull_coords")
    # If it already exists, then the pressure coordinates are converted already
    # TODO: Deal with restarts...
    isdir(pfull_dir) && return nothing
    mkdir(pfull_dir)
    if isdir(_pfull_dir)
        pfull_nc_filepaths = filter!(
            filepath -> isfile(filepath) && endswith(filepath, ".nc"),
            readdir(_pfull_dir, join = true),
        )
        for pfull_nc_filepath in pfull_nc_filepaths
            write_h_indices_to_regular_grid(
                pfull_nc_filepath,
                remapper,
                pfull_dir,
                pfull_diag_handler,
                writer,
            )
        end
    end
end

function write_h_indices_to_regular_grid(
    filepath::String,
    remapper,
    output_dir,
    pfull_diag_handler,
    netcdf_writer,
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
    target_lon, target_lat = netcdf_writer.hpts

    # Find variable name by finding the variable with 3 dimensions (time x pressure level x horizontal_index)
    not_dim_names = setdiff!(keys(ds), NCDatasets.dimnames(ds))
    varidx = findfirst(
        name ->
            NCDatasets.dimnames(ds[name]) ==
            ("time", "pressure_level", "horizontal_index"),
        not_dim_names,
    )
    varname = not_dim_names[varidx]

    data = Array(ds[varname])

    # Do interpolation here!
    interpolated_data =
        convert_to_regular_grid(ds, varname, remapper, netcdf_writer)

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
        # date_attribs = attribs_to_nt(ds["date"].attrib)
        # add_dimension!(nc_hcoords, "date", dates; date_attribs...)

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
        # deflatelevel = writer.compression_level,
        attrib = ds[varname].attrib,
    )
    v[:, :, :, :] = interpolated_data

    close(ds)
    close(nc_hcoords)
    return nothing
end



function create_pfull_coords_remapper(
    pfull_diag_handler::PfullCoordsDiagnosticsHandler,
    netcdf_writer::NetCDFWriter,
)
    # TODO: Determine if this should be on CPU or GPU (not sure if it is worth transfering the data to
    # GPU...)
    # Get pressure field and any level of the pressure field
    pfull_field = pfull_diag_handler.pfull_field
    pfull_field_level = ClimaCore.to_cpu(ClimaCore.level(pfull_field, 1))

    # Create Remapper object with specified horizontal and vertical points
    space = axes(pfull_field_level)
    lons, lats = netcdf_writer.hpts
    target_hcoords =
        [Geometry.LatLongPoint(lat, lon) for lon in lons, lat in lats]
    remapper = Remapping.Remapper(space, target_hcoords)
    return (; remapper, pfull_field_level)
end

# Possible implementation of conversion to regular grid
# Not sure if I want to pass in a netcdf file or open it myself...
function convert_to_regular_grid(
    ds,
    varname,
    remapper,
    netcdf_writer::NetCDFWriter,
)
    (; remapper, pfull_field_level) = remapper
    # Note: Order is ("time", "pressure_level", "lon", "lat")
    lons, lats = netcdf_writer.hpts
    lon_length = length(lons)
    lat_length = length(lats)
    pfull_length = length(ds["pressure_level"])
    time_length = length(ds["time"])

    interpolated_data = zeros(time_length, pfull_length, lon_length, lat_length)
    data = Array(ds[varname])

    for (idx, level) in pairs(eachslice(data, dims = (1, 2)))
        values = ClimaCore.DataLayouts.array2data(
            level,
            ClimaCore.Fields.field_values(pfull_field_level),
        )
        field = ClimaCore.Fields.Field(values, axes(pfull_field_level))
        dest = @view interpolated_data[idx, :, :]
        Remapping.interpolate!(dest, remapper, field)
    end
    return interpolated_data
end

end
