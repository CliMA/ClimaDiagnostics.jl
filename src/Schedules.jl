"""
The `Schedules` module also contain a collection of predefined schedule
functions to implement the most common behaviors (e.g., run the callback every N
steps).
"""
module Schedules

import Dates

import ..seconds_to_str_short,
    ..seconds_to_str_long,
    ..time_to_date,
    ..period_to_str_short,
    ..period_to_str_long

import SciMLBase

"""
    AbstractSchedule

`AbstractSchedule`s are structs that behave like functions and are used for the purpose of
defining a schedule to be used in `ScheduledDiagnostics`. They also may contain additional
information.
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
function (schedule::DivisorSchedule)(integrator)::Bool
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
    t_last::Base.RefValue{T}

    """The interval of time needed to elapse for the next time that this function will
    return true."""
    dt::T

    """
        EveryDtSchedule(dt; t_start = zero(dt))

    True every time the current time is larger than the previous time this schedule was true + dt.
    """
    function EveryDtSchedule(dt::T; t_start::T = zero(dt)) where {T}
        new{typeof(dt)}(Ref(t_start), dt)
    end
end

"""
    EveryDtSchedule(integrator)

Returns true if `integrator.t >= last_t + dt`, where `last_t` is the last time
this function was true and `dt` is the schedule interval time.
"""
function (schedule::EveryDtSchedule)(integrator)::Bool
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


"""
    EveryCalendarDtSchedule

Returns true if `dt` has passed since the last time this schedule was true. `dt`
here is a `Dates.Period` (e.g., `Dates.Month(1)`).

!!! compat "ClimaDiagnostics 0.2.4"
    This schedule was introduced in version `0.2.4`.
"""
struct EveryCalendarDtSchedule{P <: Dates.Period, T <: AbstractFloat} <:
       AbstractSchedule
    """Last date this function returned true."""
    date_last::Base.RefValue{Dates.DateTime}

    """The `Dates.Period` needed to elapse for the next time that this function will return
    true."""
    dt::P

    """The `Dates.DateTime` used to convert from simulation time to date."""
    reference_date::Dates.DateTime

    """Simulation time at the beginning of the simulation. Typically used in restarts."""
    t_start::T

    """
        EveryCalendarDtSchedule(dt::Dates.Period,
                                reference_date::Dates.DateTime,
                                date_last::Union{Nothing, Dates.DateTime} = nothing,
                                t_start = 0)

    True every time the current date is larger than the previous date this schedule was true + dt.

    The date is computed assuming that `integrator.t` is in seconds and using `reference_date`.
    Schematically:
    ```julia
    date = reference_date + Second(t_start + integrator.t)
    ```

    When `date_last` is `nothing`, `date_last` is computed from `reference_date` and `t_start`.

    !!! compat "ClimaDiagnostics 0.2.4"
        This schedule was introduced in version `0.2.4`.
    """
    function EveryCalendarDtSchedule(
        dt::Dates.Period;
        reference_date::Union{Dates.Date, Dates.DateTime},
        date_last::Union{Dates.DateTime, Nothing} = nothing,
        t_start::AbstractFloat = 0.0,
    )
        isnothing(date_last) &&
            (date_last = time_to_date(t_start, reference_date))
        new{typeof(dt), typeof(t_start)}(
            Ref(date_last),
            dt,
            reference_date,
            t_start,
        )
    end
end

"""
    EveryCalendarDtSchedule(integrator)

Returns true if `current_date >= last_date + dt`, where `last_date` is the last time
this function was true and `dt` is the schedule interval time.

`current_date` is computed using the schedule `reference_date` and `t_start`.
See constructor for more information.
"""
function (schedule::EveryCalendarDtSchedule)(integrator)::Bool
    next_date = schedule.date_last[] + schedule.dt
    reference_date = schedule.reference_date
    current_date = time_to_date(integrator.t, reference_date)
    if current_date >= next_date
        schedule.date_last[] = current_date
        return true
    else
        return false
    end
end

"""
    short_name(schedule::EveryCalendarDtSchedule)

Short of name of the given `schedule`. Typically used in names of files/datasets.

By default, the name of this schedule is the value converted into DDd_HHh_MMm_SSs.

Note:

This assumes that units are seconds.
"""
function short_name(schedule::EveryCalendarDtSchedule)
    return period_to_str_short(schedule.dt)
end

"""
    long_name(schedule::EveryCalendarDtSchedule)

Long of name of the given `schedule`. Typically used in attributes.

This is directly the string representation of a `Dates.Period`.
"""
function long_name(schedule::EveryCalendarDtSchedule)
    return period_to_str_long(schedule.dt)
end

end
