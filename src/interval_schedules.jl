# Interval-based output schedules. A schedule can report the accumulation window `[lo, hi)` a
# reduction's samples belong to (via `should_accumulate`/`output_period_bounds`, which have
# generic fallbacks so existing schedules are unchanged). `include`d into `module Schedules`.

"""
    ScheduleInterval(lo, hi)

A left-closed, right-open window `[lo, hi)` of `Dates.DateTime`s reported by
[`CalendarIntervalSchedule`](@ref). Reductions are timestamped at `lo` and output at `hi`.
"""
struct ScheduleInterval{T}
    lo::T
    hi::T
end

Base.show(io::IO, i::ScheduleInterval) = print(io, "[$(i.lo), $(i.hi))")
Base.:(==)(a::ScheduleInterval, b::ScheduleInterval) =
    a.lo == b.lo && a.hi == b.hi

"""
    should_accumulate(schedule, integrator)

Whether the current time lies in the *currently open* active accumulation window. A reduction
folds a sample only when its compute schedule fires and this is `true`. Defaults to `true`.

!!! note "Called around the schedule's output call"
    For stateful schedules (e.g. [`CalendarIntervalSchedule`](@ref)) the open window advances
    when the schedule's output call fires, so `orchestrate_diagnostics` reads this twice: once
    *before* the output call (does the sample fold into the closing window?) and, when the
    schedule fired and the first read was `false`, once *after* it (does the boundary sample
    belong to the newly opened window?). In the latter case the sample is folded after the
    output is written and the accumulator reset.

!!! note "Called at construction"
    `DiagnosticsHandler` also calls this once at setup with a minimal `(; t)` integrator (only
    `t` is available) to decide whether to seed the `t0` sample. Custom implementations should
    rely only on `integrator.t` (or the schedule's own state), not other integrator fields.
"""
should_accumulate(schedule, integrator) = true

"""
    output_period_bounds(schedule)

The `[lo, hi)` bounds of the most recently *closed and active* window of `schedule`, or
`nothing` if the schedule reports no window (the default) or no active window has closed yet.
When non-`nothing`, the NetCDF writer stamps these exact bounds in `time_bnds`/`date_bnds`
for reductions; otherwise it reconstructs the bounds as before. This lets the writer use a
schedule's bounds without changing the `write_field!` signature. Window reporting applies to
reductions only -- instantaneous diagnostics always use reconstructed bounds.
"""
output_period_bounds(schedule) = nothing

"""
    CalendarIntervalSchedule(period; start_date, date_last = start_date, spinup = nothing, name = "")
    CalendarIntervalSchedule(period, t::ITime; spinup = nothing, name = "")

Output schedule that partitions time into calendar windows of length `period` (a
`Dates.Period`, e.g. `Dates.Month(1)`) and closes the window `[date_last, date_last + period)`
when the date reaches a `period` boundary, like [`EveryCalendarDtSchedule`](@ref). It also
reports the closed window as a `ScheduleInterval` (so the NetCDF writer writes those exact
`time_bnds`/`date_bnds`) and, if `spinup` is set, marks windows starting before `spinup`
inactive (no accumulation, no output).

Windows are exact left-closed partitions: a sample at time `t` belongs to the window with
`lo <= date(t) < hi`, so the sample at a boundary (including `t0`) is folded into the window
*starting* there. Consequently, restarting a simulation exactly at a closed boundary
reproduces an uninterrupted run exactly: the boundary sample is re-seeded into the new window
at `DiagnosticsHandler` construction, just as an uninterrupted run folds it after the output.

`start_date` converts integrator time to a date; for `ITime` integrators it must equal the
epoch (the second constructor sets `start_date = epoch(t)` and `date_last = date(t)`).

Windows follow calendar arithmetic and may have unequal lengths when `start_date` is not the
first of the period (e.g. a window opened on Jan 31 closes on Feb 28/29). The exact closed
window is reported via [`output_period_bounds`](@ref) and written to `time_bnds`/`date_bnds`
for reductions; instantaneous diagnostics ignore the window and use reconstructed bounds.
A step that overshoots a boundary (`hi < t < hi + dt`) still closes the exact calendar window;
its sample is attributed to the *next* window.

For NetCDF output, `start_date` must equal the `NetCDFWriter`'s `start_date` (the simulation
`t = 0` date): the writer converts the reported window to seconds relative to *its own*
`start_date`, so a mismatch silently shifts the `time` axis. Use one `start_date` everywhere.

The schedule advances by at most one `period` per call, so `dt` must be much smaller than
`period` (the usual case). If a single step spans two or more periods, only one window closes
per step and output lags real time.

The struct is stateful: it mutates `date_last` and the last-closed window, so do not share an
instance between diagnostics (that would interleave their window state). On restart, set
`date_last` to the last closed window boundary at or before the restart date, or use the
`ITime` form (which sets `date_last` to the restart date and is therefore only exact when the
restart is at a window boundary).

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

# Exact left-closed membership: a sample folds into the open window `[date_last,
# date_last + period)` only if its date lies inside it. A boundary sample fails this test
# before the schedule advances (it is past `hi`) and passes it afterwards (it starts the new
# window), which is how `orchestrate_diagnostics` defers it past the output.
function should_accumulate(s::CalendarIntervalSchedule, integrator)
    _window_active(s) || return false
    current_date = time_to_date(integrator.t, s.start_date)
    return s.date_last[] <= current_date < s.date_last[] + s.period
end

"""
    (schedule::CalendarIntervalSchedule)(integrator)

Return `true` and record the closed window when the date reaches a `period` boundary (unless
the window is inactive). This call advances the open window, so `should_accumulate` reads the
closing window before it and the newly opened one after it.
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

"""
    IntervalSchedule(schedule; start_date, date_last = start_date, spinup = nothing, name = "")
    IntervalSchedule(schedule, t::ITime; spinup = nothing, name = "")

Wrap any output `schedule` into an interval schedule: the wrapped schedule's firing times
become the window boundaries. If the wrapper opens its first window at `date_last` and the
wrapped schedule fires at `t1, t2, ...`, the reported windows are `[date_last, t1)`,
`[t1, t2)`, ... -- an exact left-closed partition with the same boundary-sample deferral,
`spinup` gating, and exact `time_bnds`/`date_bnds` reporting as
[`CalendarIntervalSchedule`](@ref).

Boundaries are the *realized* firing times: a step that overshoots the wrapped schedule's
nominal boundary becomes the window edge, so edges follow the stepper, not a nominal grid
(for calendar grids use `CalendarIntervalSchedule`). On restart, set `date_last` to the date
of the last firing; this is exact when the restart coincides with a firing time.

`start_date` converts integrator time to dates and, as for `CalendarIntervalSchedule`, must
equal the `NetCDFWriter`'s `start_date`.

To classify a boundary sample before the output call, the wrapper consults the wrapped
schedule from within [`should_accumulate`](@ref) (once per step, cached by time). The wrapped
schedule must therefore not be shared with or called by anything else. It is also consulted
once at `DiagnosticsHandler` construction with a minimal `(; t)` integrator, so it must only
read `integrator.t` (true for `EveryDtSchedule` and `EveryCalendarDtSchedule`, but not for
step-counting schedules like `DivisorSchedule`/`EveryStepSchedule`).

A firing at the very start of the open window closes a zero-length window, which is dropped
(nothing is written).

!!! compat "ClimaDiagnostics 0.3.7"
    This schedule was introduced in version `0.3.7`.
"""
struct IntervalSchedule{S} <: AbstractSchedule
    """Wrapped schedule whose firing times define the window boundaries."""
    schedule::S
    """Start of the currently open window; advances at each firing."""
    date_last::Base.RefValue{Dates.DateTime}
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
    """Time of the last wrapped-schedule consultation (one per step)."""
    time_seen::Base.RefValue{Any}
    """Whether the wrapped schedule fired at `time_seen` (consumed by the output call)."""
    fired::Base.RefValue{Bool}
end

function IntervalSchedule(
    schedule;
    start_date::Dates.DateTime,
    date_last::Dates.DateTime = start_date,
    spinup::Union{Nothing, Dates.DateTime} = nothing,
    name::AbstractString = "",
)
    return IntervalSchedule(
        schedule,
        Base.RefValue(date_last),
        start_date,
        spinup,
        String(name),
        Base.RefValue{Union{Nothing, ScheduleInterval{Dates.DateTime}}}(
            nothing,
        ),
        Base.RefValue{Any}(nothing),
        Base.RefValue(false),
    )
end

IntervalSchedule(
    schedule,
    t::ITime;
    spinup::Union{Nothing, Dates.DateTime} = nothing,
    name::AbstractString = "",
) = IntervalSchedule(
    schedule;
    start_date = epoch(t),
    date_last = date(t),
    spinup,
    name,
)

_window_active(s::IntervalSchedule) =
    isnothing(s.spinup) || s.date_last[] >= s.spinup

# Consult the wrapped schedule at most once per time: the first call at a new `t` advances it
# and caches whether it fired; later calls at the same `t` reuse the cached answer.
function _consult!(s::IntervalSchedule, integrator)
    isequal(s.time_seen[], integrator.t) && return nothing
    s.time_seen[] = integrator.t
    s.fired[] = s.schedule(integrator)
    return nothing
end

# A firing time is a window boundary, so the sample there belongs to the *next* window: not
# accumulated before the output call (`fired` still set), accumulated after it (consumed).
function should_accumulate(s::IntervalSchedule, integrator)
    _consult!(s, integrator)
    s.fired[] && return false
    return _window_active(s)
end

"""
    (schedule::IntervalSchedule)(integrator)

Return `true` and record the closed window `[date_last, now)` when the wrapped schedule fires
(unless the window is inactive or zero-length). This call advances the open window, so
`should_accumulate` reads the closing window before it and the newly opened one after it.
"""
function (s::IntervalSchedule)(integrator)
    _consult!(s, integrator)
    s.fired[] || return false
    s.fired[] = false
    lo = s.date_last[]
    hi = time_to_date(integrator.t, s.start_date)
    active = _window_active(s)
    s.date_last[] = hi
    if active && hi > lo
        s.last_interval[] = ScheduleInterval(lo, hi)
        return true
    end
    return false
end

output_period_bounds(s::IntervalSchedule) = s.last_interval[]

function short_name(s::IntervalSchedule)
    isempty(s.name) || return s.name
    base =
        s.schedule isa AbstractSchedule ? short_name(s.schedule) : "interval"
    isnothing(s.spinup) ? base :
    base * "_spinup" * Dates.format(s.spinup, "yyyymmdd")
end

function long_name(s::IntervalSchedule)
    isempty(s.name) || return s.name
    base = s.schedule isa AbstractSchedule ? long_name(s.schedule) : "interval"
    isnothing(s.spinup) ? base :
    base * ", spinup from " * Dates.format(s.spinup, "yyyy-mm-dd")
end
