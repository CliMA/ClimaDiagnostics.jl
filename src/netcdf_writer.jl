import Dates

import ClimaCore: Domains, Geometry, Grids, Fields, Meshes, Spaces
import ClimaCore.Remapping: Remapper, interpolate, interpolate!
import ..Schedules: EveryStepSchedule
import ClimaUtilities.TimeManager: ITime, date

import NCDatasets
import NVTX

# Defines target_coordinates, add_space_coordinates_maybe!, add_time_maybe! for a bunch of
# Spaces
include("netcdf_writer_coordinates.jl")

"""
    NetCDFWriter

A struct to remap `ClimaCore` `Fields` to rectangular grids and save the output to NetCDF
files.
"""
struct NetCDFWriter{
    T,
    TS,
    DI,
    SYNC,
    ZSM <: Union{AbstractZSamplingMethod, Nothing},
    DATE,
    HPTS,
    VPTS,
    GA <: Union{AbstractDict{String, String}, Nothing},
    TIME,
} <: AbstractWriter
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

    """Instance of a type that determines how the points along the vertical direction are
    sampled."""
    z_sampling_method::ZSM

    """Areas of memory preallocated where the interpolation output is saved. Only the root
    process uses this."""
    preallocated_output_arrays::DI

    """Callable that determines when to call NetCDF.sync. NetCDF.sync is needed to flush the
    output to disk. Usually, the NetCDF is able to determine when to write the output to
    disk, but this sometimes fails (e.g., on GPUs). In that case, it is convenient to force
    NetCDF to write to disk. When `sync_schedule = nothing`, it is up to NetCDF to manage
    when to write to disk. Alternatively, pass a `schedule`, a function that takes the
    integrator as input and returns a boolean"""
    sync_schedule::SYNC

    """Set of datasets that need to be synced. Useful when `sync_schedule` is not `nothing`."""
    unsynced_datasets::Set{NCDatasets.NCDataset}

    """Date of the beginning of the simulation (it is used to convert seconds to dates)."""
    start_date::DATE

    """The horizontal points along the long-lat dimensions or x-y dimensions."""
    hpts::HPTS

    """The vertical points along the vertical dimension."""
    vpts::VPTS

    """Global attributes in each NetCDF file"""
    global_attribs::GA

    """Initial time of the simulation"""
    init_time::TIME

    # TODO: Add option to write dates as time
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
    NetCDFWriter(space, output_dir)

Save a `ScheduledDiagnostic` to a NetCDF file inside the `output_dir` of the simulation by
performing a pointwise (non-conservative) remapping first.

Keyword arguments
==================

- `space`: `Space` where the `Fields` are defined. This is the most general space across the
           `Fields`. In general, this is a 3D space. From a 3D space, you can take slices and
            write 2D Fields, but the opposite is not true.
- `output_dir`: The base folder where the files should be saved.
- `num_points`: How many points to use along the different dimensions to interpolate the
                fields. This is a tuple of integers, typically having meaning Long-Lat-Z,
                or X-Y-Z (the details depend on the configuration being simulated).
- `z_sampling_method`: Instance of a `AbstractZSamplingMethod` that determines how points
                       on the vertical direction should be chosen. By default, the vertical
                       points are sampled on the grid levels.
- `compression_level`: How much to compress the output NetCDF file (0 is no compression, 9
                       is maximum compression).
- `sync_schedule`: Schedule that determines when to call `NetCDF.sync` (to flush the output
                   to disk). When `NetCDF.sync` is called, you can guarantee that the bits
                   are written to disk (instead of being buffered in memory). A schedule is
                   a boolean callable that takes as a single argument the `integrator`.
                   `sync_schedule` can also be set as `nothing`, in which case we let
                   handling buffered writes to disk.
- `start_date`: Date of the beginning of the simulation.
- `horizontal_pts`: A tuple of vectors of floats meaning Long-Lat or X-Y (the
  details depend on the configuration being simulated).
- `global_attribs`: Optional dictionary of global attributes to include in all NetCDF files
                    produced by this `NetCDFWriter`. These attributes are useful for storing
                    metadata such as `source`, `creation_date`, or `frequency`. Must be
                    `nothing` or a subtype of `AbstractDict{String, String}`. Default is
                    `nothing`.
- `init_time`: The time that the simulation is initialized with. The default is `0.0`. The
               initial time of the simulation does not need to be zero. For instance, when
               restarting a simulation, the initial time of the simulation is non-zero. If
               the simulation does not begin at `t = 0` and nothing is passed in, then the
               result could be wrong.
"""
function NetCDFWriter(
    space::Spaces.AbstractSpace,
    output_dir;
    num_points = default_num_points(space),
    compression_level = 9,
    sync_schedule = ClimaComms.device(space) isa ClimaComms.CUDADevice ?
                    EveryStepSchedule() : nothing,
    z_sampling_method = LevelsMethod(),
    start_date = nothing,
    horizontal_pts = nothing,
    global_attribs = nothing,
    init_time = 0.0,
)
    horizontal_space = Spaces.horizontal_space(space)
    is_horizontal_space = horizontal_space == space

    if is_horizontal_space
        hpts = target_coordinates(space, num_points)
        vpts = []
    else
        if z_sampling_method isa LevelsMethod
            # It is a little tricky to override the number of vertical points because we don't
            # know if the vertical direction is the 2nd (as in a plane) or 3rd index (as in a
            # box or sphere). To set this value, we check if we are on a plane or not

            # TODO: Get the number of dimensions directly from the space
            num_horiz_dimensions =
                Spaces.horizontal_space(space) isa
                Spaces.SpectralElementSpace1D ? 1 : 2

            num_vpts = Spaces.nlevels(space)

            # For any configuration, it is reasonable to assume that the last
            # value of `num_pts` is the number of vertical points
            last(num_points) != num_vpts &&
                @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
            num_points =
                Tuple([num_points[1:num_horiz_dimensions]..., num_vpts])
        end
        hpts, vpts = target_coordinates(space, num_points, z_sampling_method)
    end

    if !isnothing(horizontal_pts)
        length(hpts) != length(horizontal_pts) && error(
            "Expected horizontal_pts to be of length $(length(hpts)); found $(length(horizontal_pts)) instead",
        )
        hpts = collect(horizontal_pts)
        if is_horizontal_space
            num_points = ntuple(i -> length(hpts[i]), length(hpts))
        else
            num_points =
                (ntuple(i -> length(hpts[i]), length(hpts))..., length(vpts))
        end
    end

    hcoords = hcoords_from_horizontal_space(
        horizontal_space,
        Meshes.domain(Spaces.topology(horizontal_space)),
        hpts,
    )
    zcoords = Geometry.ZPoint.(vpts)
    remapper = Remapper(space, hcoords, zcoords)
    comms_ctx = ClimaComms.context(space)

    if is_horizontal_space
        interpolated_physical_z = []
    else
        coords_z = Fields.coordinate_field(space).z
        maybe_move_to_cpu =
            ClimaComms.device(coords_z) isa ClimaComms.CUDADevice &&
            ClimaComms.iamroot(comms_ctx) ? Array : identity

        interpolated_physical_z =
            maybe_move_to_cpu(interpolate(remapper, coords_z))
    end

    preallocated_arrays =
        ClimaComms.iamroot(comms_ctx) ?
        Dict{ScheduledDiagnostic, ClimaComms.array_type(space)}() :
        Dict{ScheduledDiagnostic, Nothing}()

    unsynced_datasets = Set{NCDatasets.NCDataset}()

    return NetCDFWriter{
        typeof(num_points),
        typeof(interpolated_physical_z),
        typeof(preallocated_arrays),
        typeof(sync_schedule),
        typeof(z_sampling_method),
        typeof(start_date),
        typeof(hpts),
        typeof(vpts),
        typeof(global_attribs),
        typeof(init_time),
    }(
        output_dir,
        Dict{String, Remapper}(),
        num_points,
        compression_level,
        interpolated_physical_z,
        Dict{String, NCDatasets.NCDataset}(),
        z_sampling_method,
        preallocated_arrays,
        sync_schedule,
        unsynced_datasets,
        start_date,
        hpts,
        vpts,
        global_attribs,
        init_time,
    )
end

function NetCDFWriter(
    space::Spaces.Spaces.FiniteDifferenceSpace,
    output_dir;
    num_points = default_num_points(space),
    compression_level = 9,
    sync_schedule = ClimaComms.device(space) isa ClimaComms.CUDADevice ?
                    EveryStepSchedule() : nothing,
    z_sampling_method = LevelsMethod(),
    start_date = nothing,
    global_attribs = nothing,
    init_time = 0.0,
)
    if z_sampling_method isa LevelsMethod
        num_vpts = Spaces.nlevels(ClimaCore.Spaces.center_space(space))
        num_vpts == last(num_points) ||
            @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
        num_points = (num_vpts,)
    end
    vpts = target_coordinates(space, num_points, z_sampling_method)
    target_zcoords = Geometry.ZPoint.(vpts)
    remapper = Remapper(space; target_zcoords)

    comms_ctx = ClimaComms.context(space)

    coords_z = Fields.coordinate_field(space).z
    maybe_move_to_cpu =
        ClimaComms.device(coords_z) isa ClimaComms.CUDADevice &&
        ClimaComms.iamroot(comms_ctx) ? Array : identity

    interpolated_physical_z = maybe_move_to_cpu(interpolate(remapper, coords_z))

    preallocated_arrays =
        ClimaComms.iamroot(comms_ctx) ?
        Dict{ScheduledDiagnostic, ClimaComms.array_type(space)}() :
        Dict{ScheduledDiagnostic, Nothing}()

    unsynced_datasets = Set{NCDatasets.NCDataset}()

    return NetCDFWriter{
        typeof(num_points),
        typeof(interpolated_physical_z),
        typeof(preallocated_arrays),
        typeof(sync_schedule),
        typeof(z_sampling_method),
        typeof(start_date),
        Nothing,
        typeof(vpts),
        typeof(global_attribs),
        typeof(init_time),
    }(
        output_dir,
        Dict{String, Remapper}(),
        num_points,
        compression_level,
        interpolated_physical_z,
        Dict{String, NCDatasets.NCDataset}(),
        z_sampling_method,
        preallocated_arrays,
        sync_schedule,
        unsynced_datasets,
        start_date,
        nothing,
        vpts,
        global_attribs,
        init_time,
    )
end

function NetCDFWriter(
    space::Spaces.Spaces.PointSpace,
    output_dir;
    compression_level = 9,
    sync_schedule = ClimaComms.device(space) isa ClimaComms.CUDADevice ?
                    EveryStepSchedule() : nothing,
    start_date = nothing,
    global_attribs = nothing,
    init_time = 0.0,
    kwargs...,
)
    comms_ctx = ClimaComms.context(space)
    preallocated_arrays =
        ClimaComms.iamroot(comms_ctx) ?
        Dict{ScheduledDiagnostic, ClimaComms.array_type(space)}() :
        Dict{ScheduledDiagnostic, Nothing}()
    unsynced_datasets = Set{NCDatasets.NCDataset}()
    return NetCDFWriter{
        Nothing,
        Nothing,
        typeof(preallocated_arrays),
        typeof(sync_schedule),
        Nothing,
        typeof(start_date),
        Nothing,
        Nothing,
        typeof(global_attribs),
        typeof(init_time),
    }(
        output_dir,
        Dict{String, Remapper}(),
        nothing,
        compression_level,
        nothing,
        Dict{String, NCDatasets.NCDataset}(),
        nothing,
        preallocated_arrays,
        sync_schedule,
        unsynced_datasets,
        start_date,
        nothing,
        nothing,
        global_attribs,
        init_time,
    )
end
"""
    interpolate_field!(writer::NetCDFWriter, field, diagnostic, u, p, t)

Perform interpolation of `field` and save output in preallocated areas of `writer`.
"""
NVTX.@annotate function interpolate_field!(
    writer::NetCDFWriter,
    field,
    diagnostic,
    u,
    p,
    t,
)

    var = diagnostic.variable

    space = axes(field)

    has_horizontal_space = !(space isa Spaces.FiniteDifferenceSpace)

    if has_horizontal_space
        horizontal_space = Spaces.horizontal_space(space)

        # We have to deal with two cases: when we have an horizontal slice (e.g., the
        # surface), and when we have a full space. We distinguish these cases by checking if
        # the given space has the horizontal_space attribute. If not, it is going to be a
        # SpectralElementSpace2D and we don't have to deal with the z coordinates.
        is_horizontal_space = horizontal_space == space
    end

    # Prepare the remapper if we don't have one for the given variable. We need one remapper
    # per variable (not one per diagnostic since all the time reductions return the same
    # type of space).

    # TODO: Expand this once we support spatial reductions.
    # TODO: More generally, this can be clean up to have less conditionals
    # depending on the type of space and use dispatch instead
    if !haskey(writer.remappers, var.short_name)

        # hpts, vpts are ranges of numbers
        # target_hcoords, target_zcoords are ranges of Geometry.Points

        target_zcoords = nothing
        target_hcoords = nothing

        if has_horizontal_space
            if is_horizontal_space
                hpts = writer.hpts
                vpts = []
            else
                hpts, vpts = writer.hpts, writer.vpts
            end

            target_hcoords = hcoords_from_horizontal_space(
                horizontal_space,
                Meshes.domain(Spaces.topology(horizontal_space)),
                hpts,
            )
        else
            vpts = writer.vpts
        end

        target_zcoords = Geometry.ZPoint.(vpts)

        writer.remappers[var.short_name] =
            Remapper(space, target_hcoords, target_zcoords)
    end

    remapper = writer.remappers[var.short_name]

    # Now we can interpolate onto the target points
    # There's an MPI call in here (to aggregate the results)
    #
    # The first time we call this, we call interpolate and allocate a new array.
    # Future calls are in-place
    if haskey(writer.preallocated_output_arrays, diagnostic)
        interpolate!(
            writer.preallocated_output_arrays[diagnostic],
            remapper,
            field,
        )
    else
        writer.preallocated_output_arrays[diagnostic] =
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
NVTX.@annotate function write_field!(
    writer::NetCDFWriter,
    field,
    diagnostic,
    u,
    p,
    t,
)
    # Only the root process has to write
    ClimaComms.iamroot(ClimaComms.context(field)) || return nothing

    space = axes(field)

    maybe_move_to_cpu =
        ClimaComms.device(field) isa ClimaComms.CUDADevice ? Array : identity
    interpolated_field =
        maybe_move_to_cpu(writer.preallocated_output_arrays[diagnostic])

    if space isa Spaces.PointSpace
        # If the space is a point space, we have to remove the singleton dimension
        interpolated_field = interpolated_field[]
    end

    FT = Spaces.undertype(space)

    output_path =
        joinpath(writer.output_dir, "$(output_short_name(diagnostic)).nc")

    if !haskey(writer.open_files, output_path)
        # Append or write a new file
        open_mode = isfile(output_path) ? "a" : "c"
        if isnothing(writer.global_attribs)
            ds = NCDatasets.Dataset(output_path, open_mode)
        else
            ds = NCDatasets.Dataset(
                output_path,
                open_mode,
                attrib = writer.global_attribs,
            )
        end
        writer.open_files[output_path] = ds
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

    dim_names = add_space_coordinates_maybe!(
        nc,
        space,
        writer.num_points,
        writer.hpts,
        writer.vpts;
        writer.z_sampling_method,
        writer.interpolated_physical_z,
    )

    start_date = get_start_date(writer, p)

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

    (v, temporal_size) = get_var_and_t_index!(
        nc,
        FT,
        writer,
        interpolated_field,
        diagnostic,
        dim_names,
        start_date,
    )

    # We need to write to the next position after what we read from the data (or the first
    # position ever if we are writing the file for the first time)
    time_index = temporal_size + 1

    # Time handling for reduced vs instantaneous diagnostics:
    # - For reduced diagnostics: store time as the START of the reduction
    #   period, with time_bnds showing [start, end] of the period.
    # - For instantaneous diagnostics: store time as the current time, with
    #   time_bnds showing [previous_time, current_time].
    isa_time_reduction = !isnothing(diagnostic.reduction_time_func)
    (; init_time) = writer
    append_temporal_values!(
        nc,
        isa_time_reduction,
        t,
        start_date,
        init_time,
        time_index,
    )

    colons = ntuple(_ -> Colon(), length(dim_names))
    v[time_index, colons...] = interpolated_field

    # Add file to list of files that might need manual sync
    push!(writer.unsynced_datasets, nc)

    return nothing
end

"""
    get_start_date(writer::NetCDFWriter, p)

Get the start date for the simulation.
"""
function get_start_date(writer::NetCDFWriter, p)
    if isnothing(writer.start_date) && hasproperty(p, :start_date)
        return getproperty(p, :start_date)
    end
    return writer.start_date
end

"""
    get_var_and_t_index!(
        nc,
        FT,
        writer::NetCDFWriter,
        interpolated_field,
        diagnostic,
        dim_names,
        start_date,
    )

Return the `var`iable stored in `nc` and the current temporal size, or if the
`var`iable does not exist, then initialize and store the variable in `nc` with
attributes about the `var`iable, and return the new NetCDF variable and 0 as the
temporal size.
"""
function get_var_and_t_index!(
    nc,
    FT,
    writer::NetCDFWriter,
    interpolated_field,
    diagnostic,
    dim_names,
    start_date,
)
    var = diagnostic.variable
    if haskey(nc, "$(var.short_name)")
        # We already have something in the file
        v = nc["$(var.short_name)"]
        temporal_size, spatial_size... = size(v)
        interpolated_size = size(interpolated_field)
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
    return (v, temporal_size)
end

"""
    append_temporal_values!(
        nc,
        isa_time_reduction,
        t,
        start_date,
        init_time,
        time_index,
    )

Append temporal data to the temporal dimensions of `nc`.

The temporal dimensions are `time`, `time_bnds`, `date`, and `date_bnds`.

Time handling for reduced vs instantaneous diagnostics:
- For reduced diagnostics: store time as the START of the reduction period, with
  time_bnds showing [start, end] of the period.
- For instantaneous diagnostics: store time as the current time, with time_bnds
  showing [previous_time, current_time].
"""
function append_temporal_values!(
    nc,
    isa_time_reduction,
    t,
    start_date,
    init_time,
    time_index,
)
    append_time_values!(nc, isa_time_reduction, time_index, t, init_time)
    append_date_values!(
        nc,
        isa_time_reduction,
        time_index,
        t,
        start_date,
        init_time,
    )
    return nothing
end

"""
    append_date_values!(
        nc,
        isa_time_reduction,
        time_index,
        t,
        start_date,
        init_time,
    )

Append date values to the `date` and `date_bnds` dimension in `nc`.
"""
function append_date_values!(
    nc,
    isa_time_reduction,
    time_index,
    t,
    start_date,
    init_time,
)
    # FIXME: We are hardcoding p.start_date !
    # FIXME: We are rounding t
    if !isnothing(start_date)
        # TODO: Use ITime here
        curr_date = start_date + Dates.Millisecond(round(1000 * float(t)))
        date_type = typeof(curr_date) # not necessarily a Dates.DateTime

        if time_index == 1
            init_date =
                init_time isa ITime ? date(init_time) :
                (start_date + Dates.Millisecond(round(1000 * float(init_time))))
        end
        if isa_time_reduction
            date_to_save =
                time_index == 1 ? init_date :
                date_type(nc["date_bnds"][2, time_index - 1])
        else
            date_to_save = curr_date
        end
        nc["date"][time_index] = date_to_save
        nc["date_bnds"][:, time_index] =
            time_index == 1 ? [init_date; curr_date] :
            [date_type(nc["date_bnds"][2, time_index - 1]); curr_date]
    end
end

"""
    append_time_values!(nc, isa_time_reduction, time_index, t, init_time)

Append time values to the `time` and `time_bnds` dimension in `nc`.
"""
function append_time_values!(nc, isa_time_reduction, time_index, t, init_time)
    # TODO: Use ITime here
    if isa_time_reduction
        # For reductions, timestamp at the start of the reduction period.
        # Assume t=0 for the first write or use the end of the previous period.
        time_to_save =
            time_index == 1 ? float(init_time) :
            nc["time_bnds"][2, time_index - 1]
    else
        # For instantaneous diagnostics, use current time
        time_to_save = float(t)
    end
    nc["time"][time_index] = time_to_save
    nc["time_bnds"][:, time_index] =
        time_index == 1 ? [float(init_time); float(t)] :
        [nc["time_bnds"][2, time_index - 1]; float(t)]
    return nothing
end

"""
    sync(writer::NetCDFWriter)

Call `NCDatasets.sync` on all the files in the `writer.unsynced_datasets` list.
`NCDatasets.sync` ensures that the values are written to file.
"""
function sync(writer::NetCDFWriter)
    foreach(NCDatasets.sync, writer.unsynced_datasets)
    empty!(writer.unsynced_datasets)
    return nothing
end

function Base.show(io::IO, writer::NetCDFWriter)
    num_open_files = length(keys(writer.open_files))
    print(
        io,
        "NetCDFWriter, writing to $(writer.output_dir) ($num_open_files files open)",
    )
end
