import Dates

import ClimaCore: Domains, Geometry, Grids, Fields, Meshes, Spaces
import ClimaCore.Remapping: Remapper, interpolate, interpolate!
import ..Schedules: EveryStepSchedule

import NCDatasets

# Defines target_coordinates, add_space_coordinates_maybe!, add_time_maybe! for a bunch of
# Spaces
include("netcdf_writer_coordinates.jl")

"""
    NetCDFWriter

A struct to remap `ClimaCore` `Fields` to rectangular grids and save the output to NetCDF
files.
"""
struct NetCDFWriter{T, TS, DI, SYNC, ZSM <: AbstractZSamplingMethod, DATE} <:
       AbstractWriter
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
"""
function NetCDFWriter(
    space,
    output_dir;
    num_points = (180, 90, 50),
    compression_level = 9,
    sync_schedule = ClimaComms.device(space) isa ClimaComms.CUDADevice ?
                    EveryStepSchedule() : nothing,
    z_sampling_method = LevelsMethod(),
    start_date = nothing,
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

            num_vpts = Meshes.nelements(Grids.vertical_topology(space).mesh)

            @warn "Disabling vertical interpolation, the provided number of points is ignored (using $num_vpts)"
            num_points =
                Tuple([num_points[1:num_horiz_dimensions]..., num_vpts])
        end
        hpts, vpts = target_coordinates(space, num_points, z_sampling_method)
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
    )
end

"""
    interpolate_field!(writer::NetCDFWriter, field, diagnostic, u, p, t)

Perform interpolation of `field` and save output in preallocated areas of `writer`.
"""
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
                writer.num_points,
                writer.z_sampling_method,
            )
        end

        hcoords = hcoords_from_horizontal_space(
            horizontal_space,
            Meshes.domain(Spaces.topology(horizontal_space)),
            hpts,
        )

        # When we disable vertical_interpolation, we override the vertical points with
        # the reference values for the vertical space.
        if writer.z_sampling_method isa LevelsMethod && !is_horizontal_space
            # We need Array(parent()) because we want an array of values, not a DataLayout
            # of Points
            vpts = Array(
                parent(
                    space.grid.vertical_grid.center_local_geometry.coordinates,
                ),
            )[
                :,
                1,
            ]
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
    space = axes(field)

    maybe_move_to_cpu =
        ClimaComms.device(field) isa ClimaComms.CUDADevice ? Array : identity
    interpolated_field =
        maybe_move_to_cpu(writer.preallocated_output_arrays[diagnostic])

    if islatlonbox(
        Meshes.domain(Spaces.topology(Spaces.horizontal_space(space))),
    )
        # ClimaCore works with LatLong points, but we want to have longitude
        # first in the output, so we have to flip things
        perm = collect(1:length(size(interpolated_field)))
        perm[1:2] .= (2, 1)
        interpolated_field = permutedims(interpolated_field, perm)
    end

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
    add_time_maybe!(
        nc,
        FT;
        units = "s",
        axis = "T",
        standard_name = "time",
        long_name = "Time",
    )

    dim_names = add_space_coordinates_maybe!(
        nc,
        space,
        writer.num_points;
        writer.z_sampling_method,
        writer.interpolated_physical_z,
    )

    start_date = nothing
    if isnothing(writer.start_date)
        start_date = get(p, :start_date, nothing)
    else
        start_date = writer.start_date
    end

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
        if !isnothing(start_date) && !haskey(v.attrib, "start_date")
            v.attrib["start_date"] = string(start_date)::String
        end
        temporal_size = 0
    end

    # We need to write to the next position after what we read from the data (or the first
    # position ever if we are writing the file for the first time)
    time_index = temporal_size + 1

    nc["time"][time_index] = t

    # FIXME: We are hardcoding p.start_date !
    # FIXME: We are rounding t
    if !isnothing(start_date)
        nc["date"][time_index] = string(start_date + Dates.Millisecond(1000t))
    end

    # TODO: It would be nice to find a cleaner way to do this
    if length(dim_names) == 3
        v[time_index, :, :, :] = interpolated_field
    elseif length(dim_names) == 2
        v[time_index, :, :] = interpolated_field
    elseif length(dim_names) == 1
        v[time_index, :] = interpolated_field
    end

    # Add file to list of files that might need manual sync
    push!(writer.unsynced_datasets, nc)

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
