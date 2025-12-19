import ClimaInterpolations

import ClimaCore: Fields

# Due to the design of DiagnosticsHandler and lack of support for pressure
# coordinates in ClimaCore, it is difficult to implement a conversion to
# pressure coordinates in DiagnosticsHandler. A more reasonable solution is to
# implement the functionality of converting to pressure coordinates in
# ClimaCore, but this requires a hand written kernel or an interface to simplify
# writing the kernel. As of now, these are not realistic options, so we write
# a diagnostics handler that specifically handle writing in pressure coordinates.
"""
    PfullCoordsDiagnosticsHandler

A struct that contains the scheduled diagnostics, ancillary data, and areas of
memory needed to store and accumulate results.

This struct differs from [`DiagnosticsHandler`](@ref) as it write the
diagnostics in pressure coordinates.
"""
struct PfullCoordsDiagnosticsHandler{
    SD,
    V <: Vector{Int},
    STORAGE,
    ACC <: Dict,
    COUNT,
    MIXIN <: AbstractHandlerMixin,
} <: AbstractDiagnosticsHandler
    """An iterable with the `ScheduledDiagnostic`s that are scheduled."""
    scheduled_diagnostics::SD

    """A Vector containing keys to index into `scheduled_diagnostics`."""
    scheduled_diagnostics_keys::V

    """Container holding a potentially pre-allocated area of memory where to
    save the newly computed results. The element type is a two dimensional
    CuArray."""
    storage::STORAGE

    """Container holding a potentially pre-allocated area of memory where to
    accumulate results. The element type is a two dimensional CuArray."""
    accumulators::ACC

    """Container holding a counter that tracks how many times the given
    diagnostics was computed from the last time it was output to disk."""
    counters::COUNT

    mixin::MIXIN
end

"""
    PfullCoordsDiagnosticsHandler(
        scheduled_diagnostics,
        Y,
        p,
        t;
        dt = nothing,
    )

An object to instantiate and manage storage spaces for `ScheduledDiagnostics`.
Unlike [`DiagnosticsHandler`](@ref), the diagnostics will be written in pressure
coordinates.

Similar to `DiagnosticsHandler`, the `PfullCoordsDiagnosticsHandler` calls
`compute!(nothing, Y, p, t)` for each diagnostic, or `compute(Y, p, t)`,
whichever is available. The result is used to allocate the areas of memory for
storage and accumulation. For diagnostics without reduction,
`write_field!(output_writer, result, diagnostic, Y, p, t)` is called too.

Unlike `DiagnosticsHandler`, `pfull_compute!` is needed to compute the pressure
field at every necessary step for saving the diagnostics in pressure
coordinates.

Note: initializing a `PfullCoordsDiagnosticsHandler` can be expensive.

Keyword arguments
===================

`dt`, if passed, is used for error checking, to ensure that the diagnostics
defined as given a given period are integer multiples of the timestep.

`pfull_levels`, if passed, are the pressure levels to interpolate to. The
pressure levels must be sorted with either forward or reverse ordering.

!!! warning "File size" # TODO: Verify if this is the case
    Due to the current implementation of `PfullCoordsDiagnosticsHandler`, the
    file size of the resulting NetCDF files can be very large.
"""
function PfullCoordsDiagnosticsHandler(
    scheduled_diagnostics,
    Y,
    p,
    t;
    dt = nothing,
)
    # TODO: Add some error handling for the space in the scheduled diagnostics
    # This can be done by checking each compute_fields
    # TODO: Make sure the reshape work for cases when it is not VIJH whatever
    # This is going to involve refactoring the reshape function


    # TODO: I would like all these checks to be moved to PfullCoordsMixin...

    # Check output_writers are NetCDFWriter
    all(
        diag.output_writer isa NetCDFWriter for diag in scheduled_diagnostics
    ) || error(
        "PfullCoordsDiagnosticsHandler only supports diagnostics using NetCDFWriter",
    )


    # TODO: Should use a singleton pattern here. Make it, so there is only
    # one PfullCoordsStyle (although really only one pfull_field and compute
    # functions)
    # Or make it work with different PfullCoordsStyle objects (for different
    # netcdf writers), since only the pfull_field and compute functions need to
    # be the same or be okay with that inefficiency
    # Might need to implement a cache for that then?

    # Check all scheduled_diagnostics have the right coordinates style
    all(
        diag.output_writer.coordinates_style isa Writers.PfullCoordsStyle for
        diag in scheduled_diagnostics
    ) || error("NetCDFWriters must use PfullCoordsStyle")

    # TODO: Can be simplified maybe (pfull field is already in coordinates_style, but only
    # one of them exist though)
    first_diag = first(scheduled_diagnostics)
    pfull_field = first_diag.output_writer.coordinates_style.pressure_field
    pfull_compute! = first_diag.output_writer.coordinates_style.pfull_compute!

    all(
        diag.output_writer.coordinates_style.pressure_field === pfull_field for
        diag in scheduled_diagnostics
    ) || error(
        "There are multiple copies of the pressure field. Only initialize one PfullCoordsStyle",
    )


    # For diagnostics that perform reductions, the storage is used for the values computed
    # at each call. Reductions also save the accumulated value in accumulators.
    storage = []
    # Not all diagnostics need an accumulator, so we put them in a dictionary
    # key-ed over the diagnostic index
    accumulators = Dict{Int, Any}()
    counters = Int[]
    scheduled_diagnostics_keys = Int[]
    # TODO: If I want to make this uniform, then this can be a dict instead of a vector or store nothing in the vector to make sure the indices work out
    compute_fields = []

    # NOTE: unique requires isequal and hash to both be implemented. We don't really want to
    # do that (implement the hash). So, we roll our own `unique`. This is O(N^2) but it is run
    # only once, so it should be fine.
    seen = []
    unique_scheduled_diagnostics = []
    for x in scheduled_diagnostics
        if all(x != sd for sd in seen)
            push!(seen, x)
            push!(unique_scheduled_diagnostics, x)
        end
    end

    if length(unique_scheduled_diagnostics) != length(scheduled_diagnostics)
        @warn "Given list of diagnostics contains duplicates, removing them"
    end

    FT = eltype(pfull_field)

    typeofarray = ClimaComms.array_type(pfull_field)

    _check_dt_schedules(dt, unique_scheduled_diagnostics)

    for (i, diag) in enumerate(unique_scheduled_diagnostics)
        push!(scheduled_diagnostics_keys, i)

        out_field = compute_field(diag, Y, p, t)

        push!(compute_fields, copy(out_field))

        pfull_levels = diag.output_writer.coordinates_style.pressure_levels
        pfull_array = typeofarray{FT}(
            zeros(length(pfull_levels), Spaces.ncolumns(axes(pfull_field))),
        )
        push!(storage, copy(pfull_array))
        push!(counters, 1)
    end

    # TODO: I don't think this work, because of how the accumulation is done :(

    # TODO: Assume to be on the same space which I think is reasonable?

    perm_matrix = typeofarray{Int32}(
        zeros(Spaces.nlevels(pfull_field), Spaces.ncolumns(axes(pfull_field))),
    )

    pfull_array = sort_pressure_columns!(pfull_field, perm_matrix)
    for diag_index in 1:length(scheduled_diagnostics)
        diag = scheduled_diagnostics[diag_index]
        interpolate_field_to_pfull_coords!(
            storage[diag_index],
            compute_fields[diag_index],
            diag.output_writer.coordinates_style.pressure_coords,
            perm_matrix,
            pfull_array,
        )
    end

    for (i, diag) in enumerate(unique_scheduled_diagnostics)
        isa_time_reduction = !isnothing(diag.reduction_time_func)
        # If it is not a reduction, call the output writer as well
        if !isa_time_reduction
            interpolate_field!(diag.output_writer, storage[i], diag, Y, p, t)
            write_field!(diag.output_writer, storage[i], diag, Y, p, t)
        else
            # Add to the accumulator

            # We use similar + .= instead of copy because CUDA 5.2 does not supported nested
            # wrappers with view(reshape(view)) objects. See discussion in
            # https://github.com/CliMA/ClimaAtmos.jl/pull/2579 and
            # https://github.com/JuliaGPU/Adapt.jl/issues/21
            accumulators[i] = similar(storage[i])
            accumulators[i] .= storage[i]
        end
    end

    compute_fields = value_types(compute_fields)[compute_fields...]
    storage = value_types(storage)[storage...]
    accumulators = Dict{Int, value_types(accumulators)}(accumulators...)

    pfull_mixin =
        PfullMixin(compute_fields, pfull_compute!, pfull_field, perm_matrix)

    return PfullCoordsDiagnosticsHandler(
        unique_scheduled_diagnostics,
        scheduled_diagnostics_keys,
        storage,
        accumulators,
        counters,
        pfull_mixin,
    )
end

"""
    sort_pressure_columns!(pfull_field, permutation_matrix)

Sort the columns of `pfull_field`, store in permutation matrix in
`permutation_matrix`, and return `pfull_field` as an array of pressures to be
used in interpolation.
"""
function sort_pressure_columns!(pfull_field::Fields.Field, permutation_matrix)
    # First axis is along the z dimension and the second axis is an enumeration of the
    # columns
    reshape_to_cols(f) =
        reshape(parent(f), Spaces.nlevels(axes(f)), Spaces.ncolumns(axes(f)))
    pfull_array = reshape_to_cols(pfull_field)
    pressure_sortperm = sortperm!(permutation_matrix, pfull_array, dims = 1)
    pfull_array .= pfull_array[pressure_sortperm]
    return pfull_array
end

"""
    interpolate_field_to_pfull_coords!(
        arrays,
        fields,
        pfull_field,
        pressure_coordinates,
        permutation_matrix,
        pfull_array,
    )

Interpolate ClimaCore `fields` to pressure coordinates.

The resulting type is a two dimensional array, where the first axis is the
pressures and the second axis is an enumeration of the columns.
"""
function interpolate_field_to_pfull_coords!(
    array,
    field,
    pressure_coordinates,
    permutation_matrix,
    pfull_array,
)
    reshape_to_cols(f) =
        reshape(parent(f), Spaces.nlevels(axes(f)), Spaces.ncolumns(axes(f)))
    field = reshape_to_cols(field)
    field .= field[permutation_matrix]
    # TODO: Check what happen if not all the values are unique...
    ClimaInterpolations.Interpolation1D.interpolate1d!(
        array,
        pfull_array,
        pressure_coordinates,
        field,
        ClimaInterpolations.Interpolation1D.Linear(),
        ClimaInterpolations.Interpolation1D.Flat(),
    )
    return nothing
end

function compute_fields!(
    diagnostic_handler::AbstractDiagnosticsHandler,
    integrator,
    active_compute,
    ::PfullMixin,
)
    (; scheduled_diagnostics, scheduled_diagnostics_keys, mixin) =
        diagnostic_handler
    (; compute_fields, pfull_compute!, pfull_field, perm_matrix) = mixin
    for diag_index in scheduled_diagnostics_keys
        active_compute[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        diagnostic_handler.counters[diag_index] += 1
        # TODO: To make this the same as the other one in clima_diagnostics, I can dispatch off of the coordinate style in diag
        # TODO: Might be worth it to create a getter for
        # coordinates_style that returns no coordinates style for most of them
        # TODO: Dispatch for compute_field
        compute_field!(
            compute_fields[diag_index],
            diag,
            integrator.u,
            integrator.p,
            integrator.t,
        )
    end

    # TODO: This need to be a new function...
    #########################
    # TODO: This can probably be written as
    # prereduction_step(diagnostics_handler, integrator, active_compute)
    # Move all the relevant fields to pressure coordinates and store in storage
    if any(active_compute)
        ### Depend on the design this can be done here or before the compute_field
        pfull_compute!(pfull_field, integrator.u, integrator.p, integrator.t)
        pfull_array = sort_pressure_columns!(pfull_field, perm_matrix)
        ###
        # TODO: I think this need to be done per diagnostic if I want to reuse the exactly same code
        for diag_index in 1:length(scheduled_diagnostics)
            active_compute[diag_index] || continue
            diag = scheduled_diagnostics[diag_index]
            interpolate_field_to_pfull_coords!(
                diagnostic_handler.storage[diag_index],
                compute_fields[diag_index],
                diag.output_writer.coordinates_style.pressure_coords,
                perm_matrix,
                pfull_array,
            )
        end
    end
    return nothing
end

"""
    IntegratorWithPfullCoordsDiagnostics(
        integrator,
        scheduled_diagnostics;
        state_name = :u,
        cache_name = :p,
    )

Return a new `integrator` with diagnostics defined by `scheduled_diagnostics`.

`IntegratorWithDiagnostics` is conceptually similar to defining a `DiagnosticsHandler`,
constructing its associated `DiagnosticsCallback`, and adding such callback to a given
integrator.

The new integrator is identical to the previous one with the only difference that it has a
new callback called after all the other callbacks to accumulate/output diagnostics.

`IntegratorWithDiagnostics` ensures that the diagnostic callbacks are initialized and called
after everything else is initialized and computed.

`IntegratorWithDiagnostics` assumes that the state is `integrator.u` and the cache is
`integrator.p`. This behavior can be customized by passing the `state_name` and `cache_name`
keyword arguments.
"""
function IntegratorWithPfullCoordsDiagnostics(
    integrator,
    scheduled_diagnostics;
    state_name = :u,
    cache_name = :p,
)
    diagnostics_handler = PfullCoordsDiagnosticsHandler(
        scheduled_diagnostics,
        getproperty(integrator, state_name),
        getproperty(integrator, cache_name),
        integrator.t;
        integrator.dt,
    )

    diagnostics_callback = DiagnosticsCallback(diagnostics_handler)

    continuous_callbacks = integrator.callback.continuous_callbacks
    discrete_callbacks =
        (integrator.callback.discrete_callbacks..., diagnostics_callback)
    callback = SciMLBase.CallbackSet(continuous_callbacks, discrete_callbacks)

    Accessors.@reset integrator.callback = callback

    return integrator
end

#! format: off
"""
    era5_pressure_levels()

Return the pressure levels used by ERA5 whose units are Pa.
"""
era5_pressure_levels() = return 100.0 .* [
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000,
]
#! format: on
