module Callbacks

import ..seconds_to_str_short, ..seconds_to_str_long

import SciMLBase

"""
    Callback

A lightweight struct that contains two functions, `callback_func`, the function that has to
be called, and `schedule_func`, a boolean function that determines if it is time to call
`func` or not.
"""
struct Callback{FUNC <: Function, SCHEDULE}
    """Function to be called. It has to take one argument, the integrator."""
    callback_func::FUNC

    """Boolean function (or, more often, a callable struct, e.g., an `AbstractSchedule`) that
    determines whether `callback_func` should be called or not. It has to take one argument,
    the integrator. Most typically, only `integrator.t` or `integrator.step` are used."""
    schedule_func::SCHEDULE
end

"""
    orchestrate_diagnostics(integrator, callbacks)

Loop over all the `callbacks`, for each, check it the condition to trigger the call is met
(by calling `callback.schedule_func(integrator)`). If yes, call
`callback.callback_func(integrator)`.

`callbacks` has to be a container of `Callback`s.
"""
function orchestrate_callbacks(integrator, callbacks)
    for callback in callbacks
        if callback.schedule_func(integrator)
            callback.callback_func(integrator)
        end
    end
    return nothing
end

"""
    CallbackOrchestrator(callbacks)

Return a `SciMLBase.Callback` that executes the diagnostic callbacks specified by
`diagnostics_callbacks` according to their schedules.

`callbacks` has to be a container of `Callback`s.
"""
function CallbackOrchestrator(callbacks)
    sciml_callback(integrator) = orchestrate_callbacks(integrator, callbacks)

    # SciMLBase.DiscreteCallback checks if the given condition is true at the end of each
    # step. So, we set a condition that is always true, the callback is called at the end of
    # every step. This callback runs `orchestrate_callbacks`, which manages which
    # diagnostics functions to call
    condition = (_, _, _) -> true

    return SciMLBase.DiscreteCallback(condition, sciml_callback)
end

#############
# Schedules #
#############
"""
    AbstractSchedule

`AbstractSchedule`s are structs that behave like functions and are used for the purpose of
defining a schedule to be used in `Callback`. They also may contain additional information.
"""
abstract type AbstractSchedule end

"""
    short_name(schedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.
"""
function short_name end

"""
    long_name(schedule)

Long of name of the given `schedule`. Typically used in attributes.
"""
function long_name end

function Base.show(io::IO, schedule::AbstractSchedule)
    # This function is used in names of files/datasets
    print(io, short_name(schedule))
end

"""
    DivisorSchedule

True when the iteration number is evenly divisible by a given number.

This is roughly equivalent to: "run this call back every N steps", with the difference that
no initial offset is possible.
"""
struct DivisorSchedule <: AbstractSchedule
    """Return true when the step number is divided evenly by this number (ie, step % divisor
        == 0) """
    divisor::Int
end

"""
    DivisorSchedule(integrator)

Returns true if `integrator.step` is evenly divided by the divisor.
"""
function (schedule::DivisorSchedule)(integrator)
    return rem(integrator.step, schedule.divisor) == 0
end

"""
    short_name(schedule::DivisorSchedule)

Short name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is `<divisor>it`, with `<divisor>` the value.
"""
function short_name(schedule::DivisorSchedule)
    return "$(schedule.divisor)it"
end

"""
    long_name(schedule::DivisorSchedule)

Long name of the given `schedule`. Typically used in attributes.

By default, the name of this schedule is "every <divisor> iterations" (even this
is not technically correct...).
"""
function long_name(schedule::DivisorSchedule)
    return "every $(schedule.divisor) iterations"
end

"""
    EveryStepSchedule()

Return a schedule that executes at the end of every step.
"""
function EveryStepSchedule()
    return DivisorSchedule(1)
end

"""
    EveryDtSchedule

True every time the current time is larger than the previous time this schedule was true + Dt.

Note, this function performs no checks on whether the step is aligned with `dt` or not.
"""
struct EveryDtSchedule{T} <: AbstractSchedule
    """The integrator time the last time this function returned true."""
    t_last::Ref{T}

    """The interval of time needed to elapse for the next time that this function will
    return true."""
    dt::T

    """
        EveryDtSchedule(dt; t_start = zero(dt))

    True every time the current time is larger than the previous time this schedule was true + dt.
    """
    function EveryDtSchedule(dt; t_start = zero(dt))
        new{typeof(dt)}(Ref(t_start), dt)
    end
end

"""
    EveryDtSchedule(integrator)

Returns true if `integrator.t >= last_t + dt`, where `last_t` is the last time
this function was true and `dt` is the schedule interval time.
"""
function (schedule::EveryDtSchedule)(integrator)
    next_t = schedule.t_last[] + schedule.dt
    # Dealing with floating point precision...
    if integrator.t > next_t || integrator.t ≈ next_t
        schedule.t_last[] = integrator.t
        return true
    else
        return false
    end
end

"""
    short_name(schedule::EveryDtSchedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is the value converted into DDd_HHh_MMm_SSs.

Note:

This assumes that units are seconds.
"""
function short_name(schedule::EveryDtSchedule)
    return seconds_to_str_short(schedule.dt)
end

"""
    long_name(schedule::EveryDtSchedule)

Short of name of the given `schedule`. Typically used in attributes.

Note:

This assumes that units are seconds.
"""
function long_name(schedule::EveryDtSchedule)
    return seconds_to_str_long(schedule.dt)
end

end
