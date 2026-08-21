import Dates

import ClimaCore: Domains, Geometry, Grids, Fields, Meshes, Spaces
import ClimaCore.Remapping: Remapper, interpolate, interpolate!
import ClimaCore.Remapping: AbstractRemappingMethod, SpectralElementRemapping
import ..Schedules: EveryStepSchedule
import ClimaUtilities.TimeManager: ITime, epoch

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
    HRM <: AbstractRemappingMethod,
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

    """Horizontal interpolation method for remapping. Passed to `ClimaCore.Remapping.Remapper`."""
    horizontal_method::HRM
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
- `start_date`: `Date` or `DateTime` of the beginning of the simulation (`t = 0`).
- `horizontal_pts`: A tuple of vectors of floats meaning Long-Lat or X-Y (the details depend
                    on the configuration being simulated).
- `global_attribs`: Optional dictionary of global attributes to include in all NetCDF files
                    produced by this `NetCDFWriter`. These attributes are useful for storing
                    metadata such as `source`, `creation_date`, or `frequency`. Must be
                    `nothing` or a subtype of `AbstractDict{String, String}`. Default is
                    `nothing`.
- `horizontal_method`: Horizontal interpolation for remapping to the output grid.
  `BilinearRemapping()` (default) uses bilinear on the 2-point cell (1D/2D horizontal);
  `SpectralElementRemapping()` uses Lagrange polynomial at spectral element quadrature points.
  Passed to `ClimaCore.Remapping.Remapper`.
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
    horizontal_method::AbstractRemappingMethod = SpectralElementRemapping(),
)
    if z_sampling_method isa RealPressureLevelsMethod &&
       space != pressure_space(z_sampling_method)
        error(
            "Space for NetCDFWriter does not match space in pressure interpolator. " *
            "When using RealPressureLevelsMethod, you must pass the space from " *
            "`pressure_space(z_sampling_method)` to the NetCDFWriter constructor.",
        )
    end
    horizontal_space = Spaces.horizontal_space(space)
    is_horizontal_space = horizontal_space == space

    if is_horizontal_space
        hpts = target_coordinates(space, num_points)
        vpts = []
    else
        num_hpts = num_horizontal_points(horizontal_space, num_points)
        # For any configuration, it is reasonable to assume that the last
        # value of `num_points` is the number of vertical points
        if z_sampling_method isa LevelsMethod
            # Levels are always sampled on the center space
            num_vpts = Spaces.nlevels(Spaces.center_space(space))
        elseif z_sampling_method isa RealPressureLevelsMethod
            num_vpts = Spaces.nlevels(space)
        else
            num_vpts = last(num_points)
        end
        last(num_points) != num_vpts &&
            @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
        num_points = (num_hpts..., num_vpts)
        hpts, vpts = target_coordinates(space, num_points, z_sampling_method)
    end

    if !isnothing(horizontal_pts)
        horizontal_space isa Spaces.MultiPointSpace && error(
            "horizontal_pts is not supported for spaces of multiple columns: the columns are written as is",
        )
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
        Meshes.domain(Spaces.grid(horizontal_space)),
        hpts,
    )
    if z_sampling_method isa RealPressureLevelsMethod
        zcoords = Geometry.PPoint.(vpts)
    else
        zcoords = Geometry.ZPoint.(vpts)
    end
    remapper = Remapper(space, hcoords, zcoords; horizontal_method)
    comms_ctx = ClimaComms.context(space)

    # If RealPressureLevelsMethod is used, then the space does not have any
    # information about z, since the vertical consists of `PPoint`s
    if is_horizontal_space || z_sampling_method isa RealPressureLevelsMethod
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
        typeof(horizontal_method),
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
        horizontal_method,
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
    horizontal_method::AbstractRemappingMethod = SpectralElementRemapping(),
)
    if z_sampling_method isa RealPressureLevelsMethod &&
       space != pressure_space(z_sampling_method)
        error(
            "Space for NetCDFWriter does not match space in pressure interpolator. " *
            "When using RealPressureLevelsMethod, you must pass the space from " *
            "`pressure_space(z_sampling_method)` to the NetCDFWriter constructor.",
        )
    end
    if z_sampling_method isa LevelsMethod
        num_vpts = Spaces.nlevels(ClimaCore.Spaces.center_space(space))
        num_vpts == last(num_points) ||
            @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
        num_points = (num_vpts,)
    elseif z_sampling_method isa RealPressureLevelsMethod
        num_vpts = Spaces.nlevels(space)
        num_vpts == last(num_points) ||
            @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
        num_points = (num_vpts,)
    end
    vpts = target_coordinates(space, num_points, z_sampling_method)
    if z_sampling_method isa RealPressureLevelsMethod
        target_zcoords = Geometry.PPoint.(vpts)
    else
        target_zcoords = Geometry.ZPoint.(vpts)
    end
    remapper = Remapper(space; target_zcoords, horizontal_method)

    comms_ctx = ClimaComms.context(space)

    interpolated_physical_z = []
    if !(z_sampling_method isa RealPressureLevelsMethod)
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
        Nothing,
        typeof(vpts),
        typeof(global_attribs),
        typeof(horizontal_method),
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
        horizontal_method,
    )
end

function NetCDFWriter(
    space::Union{Spaces.PointSpace, Spaces.MultiPointSpace},
    output_dir;
    compression_level = 9,
    sync_schedule = ClimaComms.device(space) isa ClimaComms.CUDADevice ?
                    EveryStepSchedule() : nothing,
    start_date = nothing,
    global_attribs = nothing,
    kwargs...,
)
    comms_ctx = ClimaComms.context(space)
    if space isa Spaces.MultiPointSpace
        num_points = (Spaces.ncolumns(space),)
        hpts = target_coordinates(space, num_points)
    else
        hpts = nothing
        num_points = nothing
    end
    preallocated_arrays =
        ClimaComms.iamroot(comms_ctx) ?
        Dict{ScheduledDiagnostic, ClimaComms.array_type(space)}() :
        Dict{ScheduledDiagnostic, Nothing}()
    unsynced_datasets = Set{NCDatasets.NCDataset}()
    horizontal_method = SpectralElementRemapping() # no horizontal remapping for point spaces
    return NetCDFWriter{
        typeof(num_points),
        Nothing,
        typeof(preallocated_arrays),
        typeof(sync_schedule),
        Nothing,
        typeof(start_date),
        typeof(hpts),
        Nothing,
        typeof(global_attribs),
        typeof(horizontal_method),
    }(
        output_dir,
        Dict{String, Remapper}(),
        num_points,
        compression_level,
        nothing,
        Dict{String, NCDatasets.NCDataset}(),
        nothing,
        preallocated_arrays,
        sync_schedule,
        unsynced_datasets,
        start_date,
        hpts,
        nothing,
        global_attribs,
        horizontal_method,
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
                Meshes.domain(Spaces.grid(horizontal_space)),
                hpts,
            )
        else
            vpts = writer.vpts
        end

        if writer.z_sampling_method isa RealPressureLevelsMethod
            target_zcoords = Geometry.PPoint.(vpts)
        else
            target_zcoords = Geometry.ZPoint.(vpts)
        end

        writer.remappers[var.short_name] = Remapper(
            space,
            target_hcoords,
            target_zcoords;
            horizontal_method = writer.horizontal_method,
        )
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
    write_field!(
        writer::NetCDFWriter,
        field::Fields.Field,
        diagnostic,
        u,
        p,
        t,
        window_start,
    )

Save the resampled `field` produced by `diagnostic` as directed by the `writer`.

Only the root process does something here.

Note: It assumes that the field is already resampled.

The target file is determined by `output_short_name(diagnostic)`. If the target file already
exists, append to it. If not, create a new file. If the file does not contain dimensions,
they are added the first time something is written.

Time handling, where `window_start` is the start of the output window that ends at `t`:
- For reduced diagnostics: timestamps are stored at `window_start`, with time_bnds showing
  [window_start, t].
- For instantaneous diagnostics: timestamps are stored at `t`, with time_bnds showing
  [window_start, t].

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
    window_start,
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
    elseif space isa Spaces.MultiPointSpace
        # vec return the data from the columns in matching order of the coordinates
        interpolated_field = vec(interpolated_field)
    end

    FT = Spaces.undertype(space)

    if writer.z_sampling_method isa RealPressureLevelsMethod
        output_path = joinpath(
            writer.output_dir,
            "$(output_short_name(diagnostic))_pressure.nc",
        )
    else
        output_path =
            joinpath(writer.output_dir, "$(output_short_name(diagnostic)).nc")
    end
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

    start_date = get_start_date(writer, p)
    t isa ITime && check_start_date(t, start_date)
    window_start isa ITime && check_start_date(window_start, start_date)

    # See section 9.1: Features and feature types of the CF conventions for the
    # featureType attribute
    if !haskey(nc.attrib, "featureType")
        if space isa Spaces.MultiPointSpace
            nc.attrib["featureType"] = "timeSeries"
        elseif space isa Spaces.MultiColumnFiniteDifferenceSpace
            nc.attrib["featureType"] = "timeSeriesProfile"
        end
    end

    # Save as associated float if t is ITime
    TT = typeof(float(t))
    # Define time coordinate
    add_time_maybe!(
        nc,
        TT;
        time_units_and_calendar(start_date)...,
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

    # Following the CF conventions, time_bnds inherits units and calendar from time
    add_time_bounds_maybe!(nc, TT)

    (v, temporal_size) = get_var_and_t_index!(
        nc,
        FT,
        writer,
        interpolated_field,
        diagnostic,
        dim_names,
        start_date,
    )

    # Attach the per-column lat/lon auxiliary coordinates according to section
    # 9.5 "Coordinates and metadata" of the CF conventions. See appendix H.2.1
    # for an example of this
    if space isa
       Union{Spaces.MultiPointSpace, Spaces.MultiColumnFiniteDifferenceSpace} &&
       !haskey(v.attrib, "coordinates")
        v.attrib["coordinates"] = "lat lon"
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
    append_time_values!(nc, isa_time_reduction, time_index, t, window_start)

    colons = ntuple(_ -> Colon(), length(dim_names))
    v[time_index, colons...] = interpolated_field

    # Add file to list of files that might need manual sync
    push!(writer.unsynced_datasets, nc)

    return nothing
end

"""
    get_start_date(writer::NetCDFWriter, p)

Get the start date for the simulation as a `DateTime`, or `nothing` if there is none.
"""
function get_start_date(writer::NetCDFWriter, p)
    start_date = writer.start_date
    if isnothing(start_date) && hasproperty(p, :start_date)
        start_date = getproperty(p, :start_date)
    end
    return isnothing(start_date) ? nothing : Dates.DateTime(start_date)
end

"""
    check_start_date(t::ITime, start_date)

Error if `t` is an `ITime` whose epoch is not `start_date`.
"""
function check_start_date(t::ITime, start_date)
    t_epoch = epoch(t)
    if !isnothing(t_epoch) && !isnothing(start_date) && t_epoch != start_date
        error(
            "start_date $start_date does not match the epoch $t_epoch of the time",
        )
    end
    return nothing
end

"""
    time_units_and_calendar(start_date)

Return the `units` attribute of the `time` coordinate and, if there is a `start_date`, its
`calendar` attribute. With a `start_date`, times are seconds since `start_date`, so they can
be decoded as dates.
"""
time_units_and_calendar(::Nothing) = (; units = "s")
function time_units_and_calendar(start_date::Dates.DateTime)
    reference_date = replace(string(start_date), "T" => " ")
    # proleptic_gregorian is the calendar of Dates.DateTime
    # In the future, if we support CFTime.jl with ITime, then this function
    # should be updated to support different calendar systems
    return (;
        units = "seconds since $reference_date",
        calendar = "proleptic_gregorian",
    )
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
        if haskey(nc, "z_physical")
            v.attrib["coordinates"] = "z_physical"
        end
        if !isnothing(start_date) && !haskey(v.attrib, "start_date")
            v.attrib["start_date"] = string(start_date)::String
        end
        temporal_size = 0
    end
    return (v, temporal_size)
end

"""
    append_time_values!(nc, isa_time_reduction, time_index, t, window_start)

Append time values to the `time` and `time_bnds` dimension in `nc`.
"""
function append_time_values!(
    nc,
    isa_time_reduction,
    time_index,
    t,
    window_start,
)
    # TODO: Use ITime here
    nc["time"][time_index] = isa_time_reduction ? float(window_start) : float(t)
    nc["time_bnds"][:, time_index] = [float(window_start); float(t)]
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
