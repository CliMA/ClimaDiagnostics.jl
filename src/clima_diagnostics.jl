
# This file contains:
#
# - The definition of what a DiagnosticVariable is. Conceptually, a DiagnosticVariable is a
#   variable we know how to compute from the state. We attach more information to it for
#   documentation and to reference to it with its short name. DiagnosticVariables can exist
#   irrespective of the existence of an actual simulation that is being run. Science
#   packages are encourage to define their set of pre-made `DiagnosticVariables`, for
#   example, ClimaAtmos comes with several diagnostics already defined (in the
#   `ALL_DIAGNOSTICS` dictionary).
#
# - A dictionary `ALL_DIAGNOSTICS` with all the diagnostics we know how to compute, keyed
#   over their short name. If you want to add more diagnostics, look at the included files.
#   You can add your own file if you want to define several new diagnostics that are
#   conceptually related. The dictionary `ALL_DIAGNOSTICS` should be considered an
#   implementation detail.
#
# - The definition of what a ScheduledDiagnostics is. Conceptually, a ScheduledDiagnostics is a
#   DiagnosticVariable we want to compute in a given simulation. For example, it could be
#   the temperature averaged over a day. We can have multiple ScheduledDiagnostics for the
#   same DiagnosticVariable (e.g., daily and monthly average temperatures).
#
#   We provide two types of ScheduledDiagnostics: ScheduledDiagnosticIterations and
#   ScheduledDiagnosticTime, with the difference being only in what domain the recurrence
#   time is defined (are we doing something at every N timesteps or every T seconds?). It is
#   much cleaner and simpler to work with ScheduledDiagnosticIterations because iterations
#   are well defined and consistent. On the other hand, working in the time domain requires
#   dealing with what happens when the timestep is not lined up with the output period.
#   Possible solutions to this problem include: uneven output, interpolation, or restricting
#   the user from picking specific combinations of timestep/output period. In the current
#   implementation, we choose the third option. So, ScheduledDiagnosticTime is provided
#   because it is the physically interesting quantity. If we know what is the timestep, we
#   can convert between the two and check if the diagnostics are well-posed in terms of the
#   relationship between the periods and the timesteps. In some sense, you can think of
#   ScheduledDiagnosticIterations as an internal representation and ScheduledDiagnosticTime
#   as the external interface.
#
# - A function to convert a list of ScheduledDiagnosticIterations into a list of
#   AtmosCallbacks. This function takes three arguments: the list of diagnostics and two
#   dictionaries that map each scheduled diagnostic to an area of memory where to save the
#   result and where to keep track of how many times the function was called (so that we
#   can compute stuff like averages).

import .Callbacks:
    Callback, CallbackOrchestrator, DivisorSchedule, EveryDtSchedule
import .Writers: write_field!, AbstractWriter

# We define all the known identities in reduction_identities.jl
include("reduction_identities.jl")

"""
    reset_accumulator!(accumulated_value, reduction_time_func)

Reset the `accumulated_value` to the identity of `reduction_time_func`.

It requires that a method for `identity_of_reduction` is defined for `reduction_time_func`.
These functions are defined in `reduction_identities.jl`.

Users can define their own by providing a method to `Diagnostics.identity_of_reduction`
"""
function reset_accumulator!(accumulated_value, reduction_time_func)
    # identity_of_reduction works by dispatching over operation.
    # The function is defined in reduction_identities.jl
    identity = identity_of_reduction(reduction_time_func)
    fill!(parent(accumulated_value), identity)
    return nothing
end

# When the reduction is nothing, do nothing
reset_accumulator!(_, reduction_time_func::Nothing) = nothing

"""
    accumulate!(accumulated_value, latest_computed_value, reduction_time_func)

Apply the binary `reduction_time_func` to `accumulated_value` and `latest_computed_value`
and store the output to `accumulated_value`.

For example, if `reduction_time_func = max`, this computes the max between the previous max
and the newly computed value and stores it to `accumulated_value` (so that it can be used in
the next iteration).
"""
function accumulate!(
    accumulated_value,
    latest_computed_value,
    reduction_time_func,
)
    accumulated_value .=
        reduction_time_func.(accumulated_value, latest_computed_value)
    return nothing
end

# When the reduction is nothing, we do not need to accumulate anything
accumulate!(_, _, reduction_time_func::Nothing) = nothing


function compute_callback!(
    integrator,
    accumulators,
    storage,
    diag,
    counters,
    compute!,
)
    compute!(storage, integrator.u, integrator.p, integrator.t)

    # accumulator[diag] is not defined for non-reductions
    diag_accumulator = get(accumulators, diag, nothing)

    accumulate!(diag_accumulator, storage, diag.reduction_time_func)
    counters[diag] += 1
    return nothing
end

function output_callback!(integrator, accumulators, storage, diag, counters)
    # Move accumulated value to storage so that we can output it (for reductions). This
    # provides a unified interface to pre_output_hook! and output, at the cost of an
    # additional copy. If this copy turns out to be too expensive, we can move the if
    # statement below.
    isnothing(diag.reduction_time_func) || (storage .= accumulators[diag])

    # Any operations we have to perform before writing to output? Here is where we would
    # divide by N to obtain an arithmetic average
    diag.pre_output_hook!(storage, counters[diag])

    # Write to disk
    write_field!(
        diag.output_writer,
        storage,
        diag,
        integrator.u,
        integrator.p,
        integrator.t,
    )

    # accumulator[diag] is not defined for non-reductions, in which case we return `nothing`
    # The dispatch in reset_accumulator! knows how to handle this
    diag_accumulator = get(accumulators, diag, nothing)

    # If we have a reduction, we have to reset the accumulator to its neutral state. (If we
    # don't have a reduction, we don't have to do anything, but this is handled by dispatch.)
    reset_accumulator!(diag_accumulator, diag.reduction_time_func)

    # Reset counter too
    counters[diag] = 0
    return nothing
end

struct DiagnosticsHandler{SD, STORAGE <: Dict, ACC <: Dict, COUNT <: Dict}
    """An iterable with the `ScheduledDiagnostic`s that are scheduled."""
    scheduled_diagnostics::SD

    """Dictionary that maps a given `ScheduledDiagnostic` to a potentially pre-allocated
    area of memory where to save the newly computed results."""
    storage::STORAGE

    """Dictionary that maps a given `ScheduledDiagnostic` to a potentially pre-allocated
    area of memory where to accumulate results."""
    accumulators::ACC

    """Dictionary that maps a given `ScheduledDiagnosticIterations` to the counter that
    tracks how many times the given diagnostics was computed from the last time it was
    output to disk."""
    counters::COUNT
end

"""

The `DiagnosticsHandler` initializes the diagnostics by calling the  `compute!`
function.

Note: initializing a `DiagnosticsHandler` can be expensive!

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
        # TODO: Use lazy broadcasted expressions here
        storage[diag] = variable.compute!(nothing, Y, p, t)
        counters[diag] = 1

        # If it is not a reduction, call the output writer as well
        if !isa_time_reduction
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

    return DiagnosticsHandler(
        scheduled_diagnostics,
        storage,
        accumulators,
        counters,
    )
end


"""
    DiagnosticsCallback(diagnostics_handler::DiagnosticsHandler)

Translate a `DiagnosticsHandler` into a SciML callback ready to be used.
"""
function DiagnosticsCallback(diagnostics_handler::DiagnosticsHandler)
    # TODO: We have two types of callbacks: to compute and accumulate diagnostics, and to
    # dump them to disk. At the moment, they all end up in the same place, but we might want
    # to keep them separate

    callback_arrays = map(diagnostics_handler.scheduled_diagnostics) do diag
        compute_callback =
            integrator -> begin
                compute_callback!(
                    integrator,
                    diagnostics_handler.accumulators,
                    diagnostics_handler.storage[diag],
                    diag,
                    diagnostics_handler.counters,
                    diag.variable.compute!,
                )
            end
        output_callback =
            integrator -> begin
                output_callback!(
                    integrator,
                    diagnostics_handler.accumulators,
                    diagnostics_handler.storage[diag],
                    diag,
                    diagnostics_handler.counters,
                )
            end
        [
            Callback(compute_callback, diag.compute_schedule_func),
            Callback(output_callback, diag.output_schedule_func),
        ]
    end

    return CallbackOrchestrator(vcat(callback_arrays...))
end
