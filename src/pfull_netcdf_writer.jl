"""
    interpolate_field!(
        writer::NetCDFWriter,
        array,
        diagnostic,
        u,
        p,
        t,
        ::PfullCoordsStyle
    )

TODO: Update this docstring
Move `array` on GPU/CPU to CPU array in `preallocated_output_arrays`.
"""
function interpolate_field!(
    writer::NetCDFWriter,
    array,
    diagnostic,
    u,
    p,
    t,
    ::PfullCoordsStyle
)
    # TODO: Rename this...
    # TODO: Maybe this function can be renamed to interpolate_field! (even though
    # it doesn't do that, since I need to mimic it if I want to use all of
    # orchestrate_diagnostics)
    preallocated_output_arrays =
        writer.coordinates_style.preallocated_output_arrays
    if !haskey(preallocated_output_arrays, diagnostic)
        preallocated_output_arrays[diagnostic] = Array(array)
    else
        copyto!(preallocated_output_arrays[diagnostic], array)
    end
    return nothing
end

"""
    write_field_in_pfull_coords!(writer::NetCDFWriter, diagnostic, u, p, t)

Save the resampled array produced by `diagnostic` as directed by the `writer`.

Only the root process does something here.

The target file is determined by `output_short_name(diagnostic)`. If the target file already
exists, append to it. If not, create a new file. If the file does not contain dimensions,
they are added the first time something is written.

Time handling:
- For reduced diagnostics: timestamps are stored at the START of the reduction period, with
  time_bnds showing [start, end] of the period. For the first write, t=0 is assumed; for
  subsequent writes, the end of the previous period is used.
- For instantaneous diagnostics: timestamps are stored at the current time, with time_bnds
  showing [previous_time, current_time].

Attributes are appended to the dataset:
- `short_name`
- `long_name`
- `units`
- `comments`
- `start_date`
"""
function write_field_in_pfull_coords!(writer::NetCDFWriter, diagnostic, u, p, t)
    output_arrays = writer.coordinates_style.preallocated_output_arrays
    # TODO: Not sure about this, but this could be stored as the field itself
    # if I passed pfull_compute! to it
    pfull_field = writer.coordinates_style.pressure_field
    # Only the root process has to write
    ClimaComms.iamroot(ClimaComms.context(pfull_field)) || return nothing

    field_as_array = output_arrays[diagnostic]

    var = diagnostic.variable
    space = axes(pfull_field)

    FT = Spaces.undertype(space)

    # TODO: Rename this, since we want avoid conflicts with names with DiagnosticsHandler
    # (e.g. saving both diagnostics with z and pressure coordinates)
    # I guess the variable name could also be stashed here too...
    _pfull_coords_dir = joinpath(writer.output_dir, "_pfull_coords")
    isdir(_pfull_coords_dir) || mkdir(_pfull_coords_dir)
    output_path =
        joinpath(_pfull_coords_dir, "$(output_short_name(diagnostic)).nc")

    if !haskey(writer.open_files, output_path)
        # Append or write a new file
        open_mode = isfile(output_path) ? "a" : "c"
        writer.open_files[output_path] =
            NCDatasets.Dataset(output_path, open_mode)
    end

    nc = writer.open_files[output_path]

    # Save as associated float if t is ITime
    TT = typeof(float(t))
    # Define time coordinate
    add_time_maybe!(
        nc,
        TT;
        units = "s",
        axis = "T",
        standard_name = "time",
        long_name = "Time",
        bounds = "time_bnds",
    )

    dim_names = add_space_coordinates_maybe!(writer, nc, field_as_array)

    start_date = nothing
    if isnothing(writer.start_date)
        if hasproperty(p, :start_date)
            start_date = getproperty(p, :start_date)
        end
    else
        start_date = writer.start_date
    end

    add_time_bounds_maybe!(
        nc,
        TT;
        comments = "time bounds for each time value",
        units = "s",
    )

    if !isnothing(start_date)
        add_date_maybe!(
            nc;
            units = "seconds since $start_date",
            bounds = "date_bnds",
        )
        add_date_bounds_maybe!(
            nc;
            comments = "date bounds for each date value",
            units = "seconds since $start_date",
        )
    end

    if haskey(nc, "$(var.short_name)")
        # We already have something in the file
        v = nc["$(var.short_name)"]
        temporal_size, spatial_size... = size(v)
        interpolated_size = size(field_as_array)
        spatial_size == interpolated_size ||
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
        if !isnothing(start_date) && !haskey(v.attrib, "start_date")
            v.attrib["start_date"] = string(start_date)::String
        end
        temporal_size = 0
    end

    # We need to write to the next position after what we read from the data (or the first
    # position ever if we are writing the file for the first time)
    time_index = temporal_size + 1

    # Time handling for reduced vs instantaneous diagnostics:
    # - For reduced diagnostics: store time as the START of the reduction
    #   period, with time_bnds showing [start, end] of the period.
    # - For instantaneous diagnostics: store time as the current time, with
    #   time_bnds showing [previous_time, current_time].
    isa_time_reduction = !isnothing(diagnostic.reduction_time_func)

    # TODO: Use ITime here
    if isa_time_reduction
        # For reductions, timestamp at the start of the reduction period.
        # Assume t=0 for the first write or use the end of the previous period.
        time_to_save =
            time_index == 1 ? zero(float(t)) :
            nc["time_bnds"][2, time_index - 1]
    else
        # For instantaneous diagnostics, use current time
        time_to_save = float(t)
    end

    nc["time"][time_index] = time_to_save
    nc["time_bnds"][:, time_index] =
        time_index == 1 ? [zero(float(t)); float(t)] :
        [nc["time_bnds"][2, time_index - 1]; float(t)]

    # FIXME: We are hardcoding p.start_date !
    # FIXME: We are rounding t
    if !isnothing(start_date)
        # TODO: Use ITime here
        curr_date = start_date + Dates.Millisecond(round(1000 * float(t)))
        date_type = typeof(curr_date) # not necessarily a Dates.DateTime

        if isa_time_reduction
            date_to_save =
                time_index == 1 ? start_date :
                date_type(nc["date_bnds"][2, time_index - 1])
        else
            date_to_save = curr_date
        end

        nc["date"][time_index] = date_to_save
        nc["date_bnds"][:, time_index] =
            time_index == 1 ? [start_date; curr_date] :
            [date_type(nc["date_bnds"][2, time_index - 1]); curr_date]
    end

    colons = ntuple(_ -> Colon(), length(dim_names))
    v[time_index, colons...] = field_as_array

    # Add file to list of files that might need manual sync
    push!(writer.unsynced_datasets, nc)

    return nothing
end

"""
    add_space_coordinates_maybe!(
        writer::NetCDFWriter,
        nc,
        array,
    )

Add coordinates for pressure levels and horizontal indices.

The horizontal index dimension is an enumeration of the columns for the
horizontal direction.

This also add longitude and latitudes as variables whose dimension is the
horizontal indices.
"""
function add_space_coordinates_maybe!(
    writer::NetCDFWriter,
    nc,
    array, # TODO: Rename array to something else (also I am not sure if I need this?)
)
    # TODO: This does not work for a single column
    pfull_levels = writer.coordinates_style.pressure_levels
    pfull_field = writer.coordinates_style.pressure_field
    # TODO: Rename pressure_levels to pressure_level
    haskey(nc, "horizontal_index") &&
        return ("pressure_level", "horizontal_index")

    FT = eltype(pfull_field)
    _, num_h_indices = size(array)

    # TODO: Add attributes for coordinates here!
    # TODO: Add these attributes
    # pressure_level:_FillValue = NaN ;
    # pressure_level:units = "hPa" ;
    # pressure_level:positive = "down" ;
    # TODO: Add this to PressureCoordinatesStyle
    add_dimension!(
        nc,
        "pressure_level",
        pfull_levels,
        long_name = "pressure",
        units = "Pa",
        stored_direction = "increasing",
        standard_name = "air_pressure",
    )
    add_dimension!(nc, "horizontal_index", collect(1:num_h_indices))


    # TODO: I don't think this is necessary, since most user are going to use
    # the preprocessed version anyway?
    lon = NCDatasets.defVar(
        nc,
        "lon",
        FT,
        ("horizontal_index",),
        deflatelevel = writer.compression_level,
    )
    lon.attrib["short_name"] = "lon"
    lon.attrib["long_name"] = "longitude"
    lon.attrib["units"] = "degrees_east"
    lon.attrib["comments"] = "Longitudes corresponding to horizontal indices"

    lat = NCDatasets.defVar(
        nc,
        "lat",
        FT,
        ("horizontal_index",),
        deflatelevel = writer.compression_level,
    )
    lat.attrib["short_name"] = "lat"
    lat.attrib["long_name"] = "latitude"
    lat.attrib["units"] = "degrees_north"
    lat.attrib["comments"] = "Latitudes corresponding to horizontal indices"

    reshape_to_cols(f) =
        reshape(parent(f), Spaces.nlevels(axes(f)), Spaces.ncolumns(axes(f)))

    move_to_cpu =
        ClimaComms.device(pfull_field) isa ClimaComms.CUDADevice ? Array :
        identity

    lat_array = move_to_cpu(
        vec(ClimaCore.level(Fields.coordinate_field(pfull_field).lat, 1)),
    )
    lon_array = move_to_cpu(
        vec(ClimaCore.level(Fields.coordinate_field(pfull_field).long, 1)),
    )

    lat[:] = lat_array
    lon[:] = lon_array

    return ("pressure_level", "horizontal_index")
end
