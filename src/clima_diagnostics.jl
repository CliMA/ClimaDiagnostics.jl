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
struct DiagnosticsHandler{
    SD <: Tuple,
    STORAGE <: Dict,
    ACC <: Dict,
    COUNT <: Dict,
}
    """An iterable with the `ScheduledDiagnostic`s that are scheduled."""
    scheduled_diagnostics::SD

    """Dictionary that maps a given `ScheduledDiagnostic` to a potentially pre-allocated
    area of memory where to save the newly computed results."""
    storage::STORAGE

    """Dictionary that maps a given `ScheduledDiagnostic` to a potentially pre-allocated
    area of memory where to accumulate results."""
    accumulators::ACC

    """Dictionary that maps a given `ScheduledDiagnostic` to the counter that tracks how
    many times the given diagnostics was computed from the last time it was output to
    disk."""
    counters::COUNT

    """Bitvectors that identify which diagnostics are active at the given step. This is here
    mostly to reduce inference allocations that would result from operations like filter."""
    active_compute::BitVector
    active_output::BitVector
    active_sync::BitVector
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
    storage = Dict()
    accumulators = Dict()
    counters = Dict()

    for diag in scheduled_diagnostics
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

        variable = diag.variable
        isa_time_reduction = !isnothing(diag.reduction_time_func)

        # The first time we call compute! we use its return value. All the subsequent times
        # (in the callbacks), we will write the result in place
        storage[diag] = variable.compute!(nothing, Y, p, t)
        counters[diag] = 1

        # If it is not a reduction, call the output writer as well
        if !isa_time_reduction
            interpolate_field!(diag.output_writer, storage[diag], diag, Y, p, t)
            write_field!(diag.output_writer, storage[diag], diag, Y, p, t)
        else
            # Add to the accumulator

            # We use similar + .= instead of copy because CUDA 5.2 does not supported nested
            # wrappers with view(reshape(view)) objects. See discussion in
            # https://github.com/CliMA/ClimaAtmos.jl/pull/2579 and
            # https://github.com/JuliaGPU/Adapt.jl/issues/21
            accumulators[diag] = similar(storage[diag])
            accumulators[diag] .= storage[diag]
        end
    end

    num_diagnostics = length(scheduled_diagnostics)
    active_compute = BitVector(ntuple(_ -> false, num_diagnostics))
    active_output = BitVector(ntuple(_ -> false, num_diagnostics))
    active_sync = BitVector(ntuple(_ -> false, num_diagnostics))

    return DiagnosticsHandler(
        Tuple(scheduled_diagnostics),
        storage,
        accumulators,
        counters,
        active_compute,
        active_output,
        active_sync,
    )
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
    scheduled_diagnostics = diagnostic_handler.scheduled_diagnostics

    active_compute = diagnostic_handler.active_compute
    active_output = diagnostic_handler.active_output
    active_sync = diagnostic_handler.active_sync

    # Compute
    for diag_index in 1:length(scheduled_diagnostics)
        active_compute[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        diag.variable.compute!(
            diagnostic_handler.storage[diag],
            integrator.u,
            integrator.p,
            integrator.t,
        )
        diagnostic_handler.counters[diag] += 1

        isa_time_reduction = !isnothing(diag.reduction_time_func)
        if isa_time_reduction
            diagnostic_handler.accumulators[diag] .=
                diag.reduction_time_func.(
                    diagnostic_handler.accumulators[diag],
                    diagnostic_handler.storage[diag],
                )
        end
    end

    # Pre-output (averages/interpolation)
    for diag_index in 1:length(scheduled_diagnostics)
        active_output[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        # Move accumulated value to storage so that we can output it (for reductions). This
        # provides a unified interface to pre_output_hook! and output, at the cost of an
        # additional copy. If this copy turns out to be too expensive, we can move the if
        # statement below.
        isnothing(diag.reduction_time_func) || (
            diagnostic_handler.storage[diag] .=
                diagnostic_handler.accumulators[diag]
        )

        # Any operations we have to perform before writing to output? Here is where we would
        # divide by N to obtain an arithmetic average
        diag.pre_output_hook!(
            diagnostic_handler.storage[diag],
            diagnostic_handler.counters[diag],
        )

        interpolate_field!(
            diag.output_writer,
            diagnostic_handler.storage[diag],
            diag,
            integrator.u,
            integrator.p,
            integrator.t,
        )
    end

    # Save to disk
    for diag_index in 1:length(scheduled_diagnostics)
        active_output[diag_index] || continue
        diag = scheduled_diagnostics[diag_index]

        write_field!(
            diag.output_writer,
            diagnostic_handler.storage[diag],
            diag,
            integrator.u,
            integrator.p,
            integrator.t,
        )
    end

    # Post-output clean-up
    for diag_index in 1:length(scheduled_diagnostics)
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
            fill!(parent(diagnostic_handler.accumulators[diag]), identity)
        end
        # Reset counter
        diagnostic_handler.counters[diag] = 0
    end

    return nothing
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
    update_diagnostic_handler_bitvectors(integrator, diagnostics_handler::DiagnosticsHandler)

Update the `active_{compute, update, sync}` bitvector in `diagnostics_handler`.

The `diagnostics_handler` contains three bitvectors that determine which actions should be
taken at the current iterations. They are preallocated mostly to avoid inference allocations
that result from operations with groups of `ScheduledDiagnostics`, but they can also be used
to determine if `orchestrate_diagnostics` should be called or not.

This function evaluates the various `schedule_func`s for all the `ScheduledDiagnostics` in
the `diagnostics_handler` and updates the bitvectors. This function should be called at
every step.
"""
function update_diagnostic_handler_bitvectors(integrator, diagnostics_handler)
    scheduled_diagnostics = diagnostics_handler.scheduled_diagnostics

    for index in 1:length(scheduled_diagnostics)
        diag = scheduled_diagnostics[index]
        diagnostics_handler.active_compute[index] =
            diag.compute_schedule_func(integrator)
        diagnostics_handler.active_output[index] =
            diag.output_schedule_func(integrator)
        diagnostics_handler.active_sync[index] = _needs_sync(diag, integrator)
    end

    return nothing
end

"""
    check_callback_condition(integrator, diagnostics_handler::DiagnosticsHandler)

Return true when `orchestrate_diagnostics` should be called.
"""
function check_callback_condition(integrator, diagnostics_handler)
    return any(diagnostics_handler.active_compute) ||
           any(diagnostics_handler.active_output) ||
           any(diagnostics_handler.active_sync)
end

"""
    DiagnosticsCallbacks(diagnostics_handler::DiagnosticsHandler)

Translate a `DiagnosticsHandler` into two SciML callbacks ready to be used.

The first updates internal counters in `diagnostics_handler` that check if the diagnostics
have to be computed. The second actually computes and outputs the diagnostics.

"""
function DiagnosticsCallbacks(diagnostics_handler::DiagnosticsHandler)

    # We use trivial condition to update the condition to run orchestrate_diagnostics
    trivial_condition = (_, _, _) -> true

    sciml_callback_update_diagnostic_handler_bitvectors(integrator) =
        update_diagnostic_handler_bitvectors(integrator, diagnostics_handler)

    orchestrate_condition =
        (_, _, integrator) ->
            check_callback_condition(integrator, diagnostics_handler)

    sciml_callback_orchestrate_diagnostics(integrator) =
        orchestrate_diagnostics(integrator, diagnostics_handler)

    continuous_callbacks = ()
    discrete_callbacks = (
        SciMLBase.DiscreteCallback(
            trivial_condition,
            sciml_callback_update_diagnostic_handler_bitvectors,
        ),
        SciMLBase.DiscreteCallback(
            orchestrate_condition,
            sciml_callback_orchestrate_diagnostics,
        ),
    )

    return SciMLBase.CallbackSet(continuous_callbacks, discrete_callbacks)
end

"""
    IntegratorWithDiagnostics(integrator,
                              scheduled_diagnostics)

Return a new `integrator` with diagnostics defined by `scheduled_diagnostics`.

`IntegratorWithDiagnostics` is conceptually similar to defining a `DiagnosticsHandler`,
constructing its associated `DiagnosticsCallbacks`, and adding such callbacks to a given
integrator.

The new integrator is identical to the previous one with the only difference that it has a
new callback called after all the other callbacks to accumulate/output diagnostics.

`IntegratorWithDiagnostics` ensures that the diagnostic callbacks are initialized and called
after everything else is initialized and computed.

`IntegratorWithDiagnostics` assumes that the state is `integrator.u` and the cache is
`integrator.p`. This behavior can be customized by passing the `state_name` and `cache_name`
keyword arguments.
"""
function IntegratorWithDiagnostics(integrator, scheduled_diagnostics)
    diagnostics_handler = DiagnosticsHandler(
        scheduled_diagnostics,
        integrator.u,
        integrator.p,
        integrator.t;
        integrator.dt,
    )
    diagnostics_callbacks =
        DiagnosticsCallbacks(diagnostics_handler).discrete_callbacks

    continuous_callbacks = integrator.callback.continuous_callbacks
    discrete_callbacks =
        (integrator.callback.discrete_callbacks..., diagnostics_callbacks...)
    callback = SciMLBase.CallbackSet(continuous_callbacks, discrete_callbacks)

    Accessors.@reset integrator.callback = callback

    return integrator
end
