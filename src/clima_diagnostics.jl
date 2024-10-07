import Accessors
import SciMLBase

import .Schedules: DivisorSchedule, EveryDtSchedule
import .Writers: interpolate_field!, write_field!, sync, AbstractWriter

# We define all the known identities in reduction_identities.jl
include("reduction_identities.jl")

"""
    DiagnosticsHandler

A struct that contains the scheduled diagnostics, ancillary data and areas of memory needed
to store and accumulate results.
"""
struct DiagnosticsHandler{SD, V <: Vector{Int}, STORAGE, ACC <: Dict, COUNT}
    """An iterable with the `ScheduledDiagnostic`s that are scheduled."""
    scheduled_diagnostics::SD

    """A Vector containing keys to index into `scheduled_diagnostics`."""
    scheduled_diagnostics_keys::V

    """Container holding a potentially pre-allocated
    area of memory where to save the newly computed results."""
    storage::STORAGE

    """Container holding a potentially pre-allocated
    area of memory where to accumulate results."""
    accumulators::ACC

    """Container holding a counter that tracks how many times the given
    diagnostics was computed from the last time it was output to disk."""
    counters::COUNT
end

"""
    value_types(
        data;
        value_map = unionall_type,
    )

Given `data`, return a type `Union{V...}` where `V` are the `Union` of all found types in
    the values of `data`.
"""
function value_types(data)
    ret_types = Set([])
    for k in eachindex(data)
        push!(ret_types, typeof(data[k]))
    end
    return Union{ret_types...}
end

"""
    DiagnosticsHandler(scheduled_diagnostics, Y, p, t; dt = nothing)

An object to instantiate and manage storage spaces for `ScheduledDiagnostics`.

The `DiagnosticsHandler` calls `compute!(nothing, Y, p, t)` for each diagnostic. The result
is used to allocate the areas of memory for storage and accumulation. For diagnostics
without reduction, `write_field!(output_writer, result, diagnostic, Y, p, t)` is called too.

Note: initializing a `DiagnosticsHandler` can be expensive.

Keyword arguments
===================

`dt`, if passed, is used for error checking, to ensure that the diagnostics defined as given
a given period are integer multiples of the timestep.
"""
function DiagnosticsHandler(scheduled_diagnostics, Y, p, t; dt = nothing)

    # For diagnostics that perform reductions, the storage is used for the values computed
    # at each call. Reductions also save the accumulated value in accumulators.
    storage = []
    # Not all diagnostics need an accumulator, so we put them in a dictionary key-ed over the diagnostic index
    accumulators = Dict{Int, Any}()
    counters = Int[]
    scheduled_diagnostics_keys = Int[]

    unique_scheduled_diagnostics = unique(scheduled_diagnostics)
    if length(unique_scheduled_diagnostics) != length(scheduled_diagnostics)
        @warn "Given list of diagnostics contains duplicates, removing them"
    end

    for (i, diag) in enumerate(unique_scheduled_diagnostics)
        if isnothing(dt)
            @warn "dt was not passed to DiagnosticsHandler. No checks will be performed on the frequency of the diagnostics"
        else
            if diag.compute_schedule_func isa EveryDtSchedule
                compute_dt = diag.compute_schedule_func.dt
                every_num_iteration = compute_dt / dt
                every_num_iteration ≈ round(every_num_iteration) || error(
                    "Compute dt ($compute_dt) for $(diag.output_short_name) is not an even multiple of the timestep ($dt)",
                )
            end
            if diag.output_schedule_func isa EveryDtSchedule
                output_dt = diag.output_schedule_func.dt
                every_num_iteration = output_dt / dt
                every_num_iteration ≈ round(every_num_iteration) || error(
                    "Output dt ($output_dt) for $(diag.output_short_name) is not an even multiple of the timestep ($dt)",
                )
            end
        end
        push!(scheduled_diagnostics_keys, i)

        variable = diag.variable
        isa_time_reduction = !isnothing(diag.reduction_time_func)

        # The first time we call compute! we use its return value. All the subsequent times
        # (in the callbacks), we will write the result in place. We call copy to acquire ownership
        # of the data in case compute! returned a reference.
        push!(storage, copy(variable.compute!(nothing, Y, p, t)))
        push!(counters, 1)

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
    storage = value_types(storage)[storage...]
    accumulators = Dict{Int, value_types(accumulators)}(accumulators...)

    return DiagnosticsHandler(
        unique_scheduled_diagnostics,
        scheduled_diagnostics_keys,
        storage,
        accumulators,
        counters,
    )
end

# Does the writer associated to `diag` need to be synced?
# It does only when it has a sync_schedule that is a callable and that
# callable returns true when called on the integrator
function _needs_sync(diag, integrator)
    hasproperty(diag.output_writer, :sync_schedule) || return false
    isnothing(diag.output_writer.sync_schedule) && return false
    return diag.output_writer.sync_schedule(integrator)
end

"""
    orchestrate_diagnostics(integrator, diagnostic_handler::DiagnosticsHandler)

Loop over all the `ScheduledDiagnostics` in `diagnostic_handler` and run compute
and output according to their schedule functions.
"""
function orchestrate_diagnostics(
    integrator,
    diagnostic_handler::DiagnosticsHandler,
)
    (; scheduled_diagnostics, scheduled_diagnostics_keys) = diagnostic_handler
    active_compute = Bool[]
    active_output = Bool[]
    active_sync = Bool[]

    for diag in scheduled_diagnostics
        push!(active_compute, diag.compute_schedule_func(integrator))
        push!(active_output, diag.output_schedule_func(integrator))
        push!(active_sync, _needs_sync(diag, integrator))
    end

    # Compute
    for diag_index in scheduled_diagnostics_keys
        active_compute[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        diag.variable.compute!(
            diagnostic_handler.storage[diag_index],
            integrator.u,
            integrator.p,
            integrator.t,
        )
        diagnostic_handler.counters[diag_index] += 1

        isa_time_reduction = !isnothing(diag.reduction_time_func)
        if isa_time_reduction
            diagnostic_handler.accumulators[diag_index] .=
                diag.reduction_time_func.(
                    diagnostic_handler.accumulators[diag_index],
                    diagnostic_handler.storage[diag_index],
                )
        end
    end

    # Pre-output (averages/interpolation)
    for diag_index in scheduled_diagnostics_keys
        active_output[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        # Move accumulated value to storage so that we can output it (for reductions). This
        # provides a unified interface to pre_output_hook! and output, at the cost of an
        # additional copy. If this copy turns out to be too expensive, we can move the if
        # statement below.
        isnothing(diag.reduction_time_func) || (
            diagnostic_handler.storage[diag_index] .=
                diagnostic_handler.accumulators[diag_index]
        )

        # Any operations we have to perform before writing to output? Here is where we would
        # divide by N to obtain an arithmetic average
        diag.pre_output_hook!(
            diagnostic_handler.storage[diag_index],
            diagnostic_handler.counters[diag_index],
        )

        interpolate_field!(
            diag.output_writer,
            diagnostic_handler.storage[diag_index],
            diag,
            integrator.u,
            integrator.p,
            integrator.t,
        )
    end

    # Save to disk
    for diag_index in scheduled_diagnostics_keys
        active_output[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        write_field!(
            diag.output_writer,
            diagnostic_handler.storage[diag_index],
            diag,
            integrator.u,
            integrator.p,
            integrator.t,
        )
    end

    # Post-output clean-up
    for diag_index in scheduled_diagnostics_keys
        diag = scheduled_diagnostics[diag_index]

        # First, maybe call sync for the writer. This might happen regardless of
        # whether the diagnostic was active or not (because diagnostics
        # typically share writers)
        active_sync[diag_index] && sync(diag.output_writer)

        active_output[diag_index] || continue

        # Reset accumulator
        isa_time_reduction = !isnothing(diag.reduction_time_func)
        if isa_time_reduction
            # identity_of_reduction works by dispatching over operation.
            # The function is defined in reduction_identities.jl
            identity = identity_of_reduction(diag.reduction_time_func)
            fill!(parent(diagnostic_handler.accumulators[diag_index]), identity)
        end
        # Reset counter
        diagnostic_handler.counters[diag_index] = 0
    end

    return nothing
end

"""
    DiagnosticsCallback(diagnostics_handler::DiagnosticsHandler)

Translate a `DiagnosticsHandler` into a SciML callback ready to be used.
"""
function DiagnosticsCallback(diagnostics_handler::DiagnosticsHandler)
    sciml_callback(integrator) =
        orchestrate_diagnostics(integrator, diagnostics_handler)

    # SciMLBase.DiscreteCallback checks if the given condition is true at the end of each
    # step. So, we set a condition that is always true, the callback is called at the end of
    # every step. This callback runs `orchestrate_callbacks`, which manages which
    # diagnostics functions to call
    condition = (_, _, _) -> true

    return SciMLBase.DiscreteCallback(condition, sciml_callback)
end

"""
    IntegratorWithDiagnostics(integrator,
                              scheduled_diagnostics;
                              state_name = :u,
                              cache_name = :p)

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
function IntegratorWithDiagnostics(
    integrator,
    scheduled_diagnostics;
    state_name = :u,
    cache_name = :p,
)
    diagnostics_handler = DiagnosticsHandler(
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
