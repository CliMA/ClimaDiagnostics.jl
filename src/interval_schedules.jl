# Interval-based output schedules. A schedule can report the accumulation window `(lo, hi]` a
# reduction's samples belong to (via `should_accumulate`/`output_period_bounds`, which have
# generic fallbacks so existing schedules are unchanged). `include`d into `module Schedules`.

"""
    ScheduleInterval(lo, hi)

A left-open, right-closed window `(lo, hi]` of `Dates.DateTime`s reported by
[`CalendarIntervalSchedule`](@ref). Reductions are timestamped at `lo` and output at `hi`.
"""
struct ScheduleInterval{T}
    lo::T
    hi::T
end

Base.show(io::IO, i::ScheduleInterval) = print(io, "($(i.lo), $(i.hi)]")
Base.:(==)(a::ScheduleInterval, b::ScheduleInterval) =
    a.lo == b.lo && a.hi == b.hi

"""
    should_accumulate(schedule, integrator)

Whether the current time lies in an active accumulation window. A reduction folds a sample
only when its compute schedule fires and this is `true`. Call it before the schedule's output
check, which may advance the window. Defaults to `true`.
"""
should_accumulate(schedule, integrator) = true

"""
    output_period_bounds(schedule)

The `(lo, hi]` bounds of the window `schedule` last closed, or `nothing` if it defines none
(then the NetCDF writer reconstructs the bounds as before). Lets the writer use a schedule's
bounds without changing the `write_field!` signature.
"""
output_period_bounds(schedule) = nothing

"""
    CalendarIntervalSchedule(period; start_date, date_last = start_date, spinup = nothing, name = "")
    CalendarIntervalSchedule(period, t::ITime; spinup = nothing, name = "")

Output schedule that partitions time into calendar windows of length `period` (a
`Dates.Period`, e.g. `Dates.Month(1)`) and closes the window `(date_last, date_last + period]`
when the date crosses a `period` boundary, like [`EveryCalendarDtSchedule`](@ref). It also
reports the closed window as a `ScheduleInterval` (so the NetCDF writer writes those exact
`time_bnds`/`date_bnds`) and, if `spinup` is set, marks windows starting before `spinup`
inactive (no accumulation, no output).

`start_date` converts integrator time to a date; for `ITime` integrators it must equal the
epoch (the second constructor sets that). The struct is stateful: do not share an instance
between diagnostics. On restart, set `date_last` to the restart date or use the `ITime` form.

!!! compat "ClimaDiagnostics 0.3.7"
    This schedule was introduced in version `0.3.7`.
"""
struct CalendarIntervalSchedule{P <: Dates.Period} <: AbstractSchedule
    """Start of the currently open window; advances as windows close."""
    date_last::Base.RefValue{Dates.DateTime}
    """Length of each window."""
    period::P
    """Date used to convert integrator time to a date."""
    start_date::Dates.DateTime
    """If not `nothing`, windows starting before this date are inactive."""
    spinup::Union{Nothing, Dates.DateTime}
    """Optional override of the auto-derived `short_name`."""
    name::String
    """Most recently closed window, read back by `output_period_bounds`."""
    last_interval::Base.RefValue{
        Union{Nothing, ScheduleInterval{Dates.DateTime}},
    }
end

function CalendarIntervalSchedule(
    period::Dates.Period;
    start_date::Dates.DateTime,
    date_last::Dates.DateTime = start_date,
    spinup::Union{Nothing, Dates.DateTime} = nothing,
    name::AbstractString = "",
)
    return CalendarIntervalSchedule(
        Base.RefValue(date_last),
        period,
        start_date,
        spinup,
        String(name),
        Base.RefValue{Union{Nothing, ScheduleInterval{Dates.DateTime}}}(
            nothing,
        ),
    )
end

CalendarIntervalSchedule(
    period::Dates.Period,
    t::ITime;
    spinup::Union{Nothing, Dates.DateTime} = nothing,
    name::AbstractString = "",
) = CalendarIntervalSchedule(
    period;
    start_date = epoch(t),
    date_last = date(t),
    spinup,
    name,
)

# The open window (starting at `date_last`) is active unless it begins before `spinup`.
_window_active(s::CalendarIntervalSchedule) =
    isnothing(s.spinup) || s.date_last[] >= s.spinup

should_accumulate(s::CalendarIntervalSchedule, integrator) = _window_active(s)

"""
    (schedule::CalendarIntervalSchedule)(integrator)

Return `true` and record the closed window when the date crosses a `period` boundary (unless
the window is inactive). `should_accumulate` reads `date_last` and so must be called first.
"""
function (s::CalendarIntervalSchedule)(integrator)
    current_date = time_to_date(integrator.t, s.start_date)
    if current_date >= s.date_last[] + s.period
        lo = s.date_last[]
        hi = lo + s.period
        active = _window_active(s)
        s.date_last[] = hi
        if active
            s.last_interval[] = ScheduleInterval(lo, hi)
            return true
        end
    end
    return false
end

output_period_bounds(s::CalendarIntervalSchedule) = s.last_interval[]

function short_name(s::CalendarIntervalSchedule)
    isempty(s.name) || return s.name
    base = period_to_str_short(s.period)
    isnothing(s.spinup) ? base :
    base * "_spinup" * Dates.format(s.spinup, "yyyymmdd")
end

function long_name(s::CalendarIntervalSchedule)
    isempty(s.name) || return s.name
    base = period_to_str_long(s.period)
    isnothing(s.spinup) ? base :
    base * ", spinup from " * Dates.format(s.spinup, "yyyy-mm-dd")
end
