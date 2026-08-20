using Test

import Dates

import ClimaTimeSteppers
import NCDatasets

import ClimaComms
@static if pkgversion(ClimaComms) >= v"0.6"
    ClimaComms.@import_required_backends
end
import ClimaCore

import ClimaDiagnostics
import ClimaDiagnostics.Schedules:
    CalendarIntervalSchedule,
    IntervalSchedule,
    EveryCalendarDtSchedule,
    EveryDtSchedule,
    should_accumulate,
    output_period_bounds,
    ScheduleInterval,
    long_name

import ClimaUtilities.TimeManager: ITime

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

# 2024 is a leap year: Jan = 31 days, Feb = 29 days, so monthly boundaries land at day 31
# (Feb 1) and day 60 (Mar 1).
@testset "CalendarIntervalSchedule (unit)" begin
    start_date = Dates.DateTime(2024, 1, 1)
    day = 86400.0

    # Step day-by-day as the orchestration does: `should_accumulate` is read before the
    # schedule call (membership in the closing window) and after it (membership in the newly
    # opened one); a sample in the latter only is a deferred boundary sample. `used` logs
    # whether each day's sample is folded into any window.
    function drive(schedule; ndays)
        fired = Tuple{Float64, Any, Bool}[]
        used = Bool[]
        for k in 0:ndays
            integrator = (; t = k * day)
            acc = should_accumulate(schedule, integrator)
            fire = schedule(integrator)
            deferred = !acc && should_accumulate(schedule, integrator)
            fire && push!(
                fired,
                (k * day, output_period_bounds(schedule), deferred),
            )
            push!(used, acc || deferred)
        end
        return fired, used
    end

    @testset "monthly windows (no spinup)" begin
        schedule = CalendarIntervalSchedule(Dates.Month(1); start_date)
        fired, used = drive(schedule; ndays = 70)

        @test length(fired) == 2
        # Boundary samples (days 31 and 60) are deferred into the window they start.
        @test fired[1] == (
            31 * day,
            ScheduleInterval(
                Dates.DateTime(2024, 1, 1),
                Dates.DateTime(2024, 2, 1),
            ),
            true,
        )
        @test fired[2] == (
            60 * day,
            ScheduleInterval(
                Dates.DateTime(2024, 2, 1),
                Dates.DateTime(2024, 3, 1),
            ),
            true,
        )
        @test all(used)                                      # exact partition: no sample dropped
        @test output_period_bounds(schedule) == fired[2][2]  # writer reads the last window
    end

    @testset "window bounds stay on the calendar under overshoot" begin
        # A step past a boundary still closes the exact calendar window (`hi` is the boundary,
        # not the step time), and the partition does not drift afterward.
        schedule = CalendarIntervalSchedule(Dates.Month(1); start_date)
        @test schedule((; t = 40 * day))   # Feb 10, past the Feb 1 boundary
        @test output_period_bounds(schedule) == ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 2, 1),
        )
        @test schedule.date_last[] == Dates.DateTime(2024, 2, 1)  # not Feb 10
        @test !schedule((; t = 50 * day))  # Feb 20: Mar 1 not yet reached
        @test schedule((; t = 60 * day))   # Mar 1
        @test output_period_bounds(schedule) == ScheduleInterval(
            Dates.DateTime(2024, 2, 1),
            Dates.DateTime(2024, 3, 1),
        )
    end

    @testset "names" begin
        @test "$(CalendarIntervalSchedule(Dates.Month(1); start_date))" == "1M"
        @test long_name(CalendarIntervalSchedule(Dates.Month(1); start_date)) ==
              "1 Month"
        spun = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            spinup = Dates.DateTime(2024, 2, 1),
        )
        @test "$(spun)" == "1M_spinup20240201"
        @test long_name(spun) == "1 Month, spinup from 2024-02-01"
        named =
            CalendarIntervalSchedule(Dates.Month(1); start_date, name = "DJF")
        @test "$(named)" == "DJF"
        @test long_name(named) == "DJF"
    end

    @testset "spinup gates accumulation and output" begin
        # spinup = Feb 1 makes [Jan 1, Feb 1) inactive; the first output window is [Feb 1, Mar 1).
        schedule = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            spinup = Dates.DateTime(2024, 2, 1),
        )
        fired, used = drive(schedule; ndays = 70)

        @test length(fired) == 1
        @test fired[1][1] == 60 * day
        @test fired[1][2] == ScheduleInterval(
            Dates.DateTime(2024, 2, 1),
            Dates.DateTime(2024, 3, 1),
        )
        # Days 0..30 fall in the inactive January window and are dropped; day 31 (= Feb 1)
        # starts the active February window and is deferred in, so accumulation runs from it.
        @test all(.!used[1:31])
        @test all(used[32:end])
    end

    @testset "ITime constructor anchors to the epoch" begin
        schedule = CalendarIntervalSchedule(
            Dates.Month(1),
            ITime(0, epoch = start_date),
        )
        @test schedule.start_date == start_date
        @test schedule.date_last[] == start_date
        @test isnothing(output_period_bounds(schedule))
    end

    @testset "boundary sample after an inactive window is deferred in" begin
        # Day 31 == Feb 1 is not in the (inactive) January window [Jan 1, Feb 1), but after
        # the schedule call — which closes January without output — it starts the active
        # February window and must be folded in.
        schedule = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            spinup = Dates.DateTime(2024, 2, 1),
        )
        @test should_accumulate(schedule, (; t = 31 * day)) == false
        @test !schedule((; t = 31 * day))  # closes inactive January, advances to Feb 1
        @test should_accumulate(schedule, (; t = 31 * day)) == true
    end

    @testset "non-first-of-month start (uneven calendar windows)" begin
        sd = Dates.DateTime(2024, 1, 31)
        s = CalendarIntervalSchedule(Dates.Month(1); start_date = sd)
        @test s((; t = 40 * day))                       # past the Feb boundary
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 1, 31),
            Dates.DateTime(2024, 2, 29),               # leap year, calendar arithmetic
        )
        @test s.date_last[] == Dates.DateTime(2024, 2, 29)
        @test !s((; t = 50 * day))
        @test s((; t = 75 * day))                       # past the Mar boundary
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 2, 29),
            Dates.DateTime(2024, 3, 29),
        )
    end

    @testset "year boundary" begin
        sd = Dates.DateTime(2024, 11, 1)
        s = CalendarIntervalSchedule(Dates.Month(1); start_date = sd)
        # Nov has 30 days, Dec has 31: boundaries at day 30 (Dec 1) and day 61 (Jan 1, 2025).
        @test s((; t = 30 * day))
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 11, 1),
            Dates.DateTime(2024, 12, 1),
        )
        @test s((; t = 61 * day))
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 12, 1),
            Dates.DateTime(2025, 1, 1),                # year rolls over cleanly
        )
    end

    @testset "restart resumes on the calendar grid" begin
        # Simulate a restart at Feb 15: date_last is the last closed boundary (Feb 1).
        s = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            date_last = Dates.DateTime(2024, 2, 1),
        )
        @test !s((; t = 50 * day))                      # Feb 20: Mar 1 not yet reached
        @test s((; t = 60 * day))                       # Mar 1
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 2, 1),
            Dates.DateTime(2024, 3, 1),
        )

        # Restarting exactly at a boundary is exact: the boundary sample belongs to the new
        # window, so the construction-time seed reproduces an uninterrupted run's deferral.
        r = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            date_last = Dates.DateTime(2024, 2, 1),
        )
        @test should_accumulate(r, (; t = 31 * day))    # Feb 1 ∈ [Feb 1, Mar 1)
    end

    @testset "weekly and daily periods" begin
        sd = Dates.DateTime(2024, 1, 1)
        w = CalendarIntervalSchedule(Dates.Week(1); start_date = sd)
        @test w((; t = 7 * day))
        @test output_period_bounds(w) == ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 1, 8),
        )

        d = CalendarIntervalSchedule(Dates.Day(10); start_date = sd)
        @test d((; t = 10 * day))
        @test output_period_bounds(d) == ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 1, 11),
        )
    end

    @testset "ScheduleInterval == and show" begin
        a = ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 2, 1),
        )
        b = ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 2, 1),
        )
        c = ScheduleInterval(
            Dates.DateTime(2024, 1, 2),
            Dates.DateTime(2024, 2, 1),
        )
        @test a == b
        @test a != c
        @test string(a) == "[2024-01-01T00:00:00, 2024-02-01T00:00:00)"
    end

    @testset "IntervalSchedule wraps EveryDtSchedule" begin
        # Firing times become the boundaries: EveryDtSchedule(30 days) fires on days 30 and
        # 60, so the windows are [Jan 1, Jan 31) and [Jan 31, Mar 1) (realized edges).
        s = IntervalSchedule(EveryDtSchedule(30 * day); start_date)
        fired, used = drive(s; ndays = 70)

        @test length(fired) == 2
        @test fired[1] == (
            30 * day,
            ScheduleInterval(
                Dates.DateTime(2024, 1, 1),
                Dates.DateTime(2024, 1, 31),
            ),
            true,
        )
        @test fired[2] == (
            60 * day,
            ScheduleInterval(
                Dates.DateTime(2024, 1, 31),
                Dates.DateTime(2024, 3, 1),
            ),
            true,
        )
        @test all(used)
    end

    @testset "IntervalSchedule consults the wrapped schedule once per time" begin
        s = IntervalSchedule(EveryDtSchedule(30 * day); start_date)
        @test should_accumulate(s, (; t = 0.0))
        @test should_accumulate(s, (; t = 0.0))        # cached, no double advance
        @test !should_accumulate(s, (; t = 30 * day))  # boundary sample
        @test s((; t = 30 * day))                      # closes [Jan 1, Jan 31)
        @test should_accumulate(s, (; t = 30 * day))   # now in the new window
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 1, 31),
        )
    end

    @testset "IntervalSchedule uses realized (overshooting) boundaries" begin
        s = IntervalSchedule(EveryDtSchedule(30 * day); start_date)
        @test !s((; t = 29 * day))
        @test s((; t = 32 * day))                      # first call at/past day 30
        @test output_period_bounds(s) == ScheduleInterval(
            Dates.DateTime(2024, 1, 1),
            Dates.DateTime(2024, 2, 2),                # Jan 1 + 32 days: the step time
        )
    end

    @testset "IntervalSchedule spinup" begin
        # spinup = Jan 31 makes [Jan 1, Jan 31) inactive; its day-30 boundary sample is
        # deferred into the active [Jan 31, Mar 1) window.
        s = IntervalSchedule(
            EveryDtSchedule(30 * day);
            start_date,
            spinup = Dates.DateTime(2024, 1, 31),
        )
        fired, used = drive(s; ndays = 70)

        @test length(fired) == 1
        @test fired[1] == (
            60 * day,
            ScheduleInterval(
                Dates.DateTime(2024, 1, 31),
                Dates.DateTime(2024, 3, 1),
            ),
            true,
        )
        @test all(.!used[1:30])
        @test all(used[31:end])
    end

    @testset "non-interval schedules use the generic fallbacks" begin
        # Guarantees the change is non-breaking: anything that does not report a window keeps
        # the default behavior (accumulate every step, no schedule-supplied bounds).
        cal = EveryCalendarDtSchedule(Dates.Month(1); start_date)
        @test isnothing(output_period_bounds(cal))
        @test should_accumulate(cal, (; t = 0.0)) == true
        plain = integrator -> true
        @test isnothing(output_period_bounds(plain))
        @test should_accumulate(plain, (; t = 0.0)) == true
    end
end

# End-to-end: a sample is folded only when the compute schedule fires AND the time is in an
# active window; a boundary sample is folded into the window it starts (after the output).
# With `value = t -> 1` and a `+` reduction (no hook), the output equals the number of folded
# samples. dt = 1 day, so boundaries land on steps 31 (Feb 1) and 60 (Mar 1).
@testset "CalendarIntervalSchedule (end-to-end)" begin
    const_start_date = Dates.DateTime(2024, 1, 1)
    day = 86400.0

    function calendar_integrator(;
        output_writer,
        output_short_name,
        spinup = nothing,
        reduction_time_func = (+),
        pre_output_hook! = nothing,
        value = t -> 1,
        schedule = nothing,
    )
        space = ColumnCenterFiniteDifferenceSpace()
        args, kwargs = create_problem(space; t0 = 0.0, tf = 70 * day, dt = day)
        function compute!(out, u, p, t)
            isnothing(out) && return u.my_var .* 0 .+ value(t)
            out .= value(t)
            return nothing
        end
        var = ClimaDiagnostics.DiagnosticVariable(;
            compute! = compute!,
            short_name = "C",
            long_name = "Constant",
        )
        diag = ClimaDiagnostics.ScheduledDiagnostic(;
            variable = var,
            output_writer,
            reduction_time_func,
            output_schedule_func = isnothing(schedule) ?
                                   CalendarIntervalSchedule(
                Dates.Month(1);
                start_date = const_start_date,
                spinup,
            ) : schedule,
            pre_output_hook!,
            output_short_name,
        )
        return ClimaDiagnostics.IntegratorWithDiagnostics(
            ClimaTimeSteppers.init(args...; kwargs...),
            [diag],
        )
    end

    @testset "sample gating (sum reduction)" begin
        # No spinup: [Jan 1, Feb 1) holds the t0 seed + days 1..30 (31 samples); [Feb 1, Mar 1)
        # holds the deferred day-31 boundary sample + days 32..59 (29).
        dw = ClimaDiagnostics.Writers.DictWriter()
        ClimaTimeSteppers.solve!(
            calendar_integrator(;
                output_writer = dw,
                output_short_name = "csum",
            ),
        )
        if ClimaComms.iamroot(context)
            out = dw["csum"]
            @test sort(collect(keys(out))) == [31 * day, 60 * day]
            @test parent(out[31 * day])[1] == 31
            @test parent(out[60 * day])[1] == 29
        end

        # spinup = Feb 1: January never accumulates, so only [Feb 1, Mar 1) is written (29 samples).
        dw_spin = ClimaDiagnostics.Writers.DictWriter()
        ClimaTimeSteppers.solve!(
            calendar_integrator(;
                output_writer = dw_spin,
                output_short_name = "csum_spin",
                spinup = Dates.DateTime(2024, 2, 1),
            ),
        )
        if ClimaComms.iamroot(context)
            out = dw_spin["csum_spin"]
            @test collect(keys(out)) == [60 * day]
            @test parent(out[60 * day])[1] == 29
        end
    end

    @testset "IntervalSchedule end-to-end (sum reduction)" begin
        # EveryDtSchedule(30 days) wrapped as intervals: [0, 30d) holds the t0 seed +
        # days 1..29 (30 samples); [30d, 60d) holds the deferred day-30 sample + days
        # 31..59 (30).
        dw = ClimaDiagnostics.Writers.DictWriter()
        ClimaTimeSteppers.solve!(
            calendar_integrator(;
                output_writer = dw,
                output_short_name = "isum",
                schedule = IntervalSchedule(
                    EveryDtSchedule(30 * day);
                    start_date = const_start_date,
                ),
            ),
        )
        if ClimaComms.iamroot(context)
            out = dw["isum"]
            @test sort(collect(keys(out))) == [30 * day, 60 * day]
            @test parent(out[30 * day])[1] == 30
            @test parent(out[60 * day])[1] == 30
        end
    end

    @testset "non-additive reduction with spinup (min)" begin
        # Each sample equals `t`; `min` over [Feb 1, Mar 1) is day 31, the deferred boundary
        # sample. A bad identity seed or a leaked pre-spinup sample (e.g. the t0 = 0 seed)
        # would drag the minimum to 0.
        dw = ClimaDiagnostics.Writers.DictWriter()
        ClimaTimeSteppers.solve!(
            calendar_integrator(;
                output_writer = dw,
                output_short_name = "cmin",
                reduction_time_func = min,
                spinup = Dates.DateTime(2024, 2, 1),
                value = t -> float(t),
            ),
        )
        if ClimaComms.iamroot(context)
            out = dw["cmin"]
            @test collect(keys(out)) == [60 * day]
            @test parent(out[60 * day])[1] == 31 * day
        end
    end

    @testset "non-additive reduction with spinup (max)" begin
        # Each sample equals `t`; the active window [Feb 1, Mar 1) folds days 31..59 (day 60
        # is deferred to March), so the `max` is day 59. A wrong identity seed (+Inf instead
        # of -Inf) or a leaked pre-spinup sample would corrupt the result.
        dw = ClimaDiagnostics.Writers.DictWriter()
        ClimaTimeSteppers.solve!(
            calendar_integrator(;
                output_writer = dw,
                output_short_name = "cmax",
                reduction_time_func = max,
                spinup = Dates.DateTime(2024, 2, 1),
                value = t -> float(t),
            ),
        )
        if ClimaComms.iamroot(context)
            out = dw["cmax"]
            @test collect(keys(out)) == [60 * day]
            @test parent(out[60 * day])[1] == 59 * day
        end
    end

    @testset "NetCDF time bounds with spinup" begin
        mktempdir() do output_dir
            output_dir = ClimaComms.bcast(context, output_dir)
            writer = ClimaDiagnostics.Writers.NetCDFWriter(
                ColumnCenterFiniteDifferenceSpace(),
                output_dir;
                num_points = (10,),
                start_date = const_start_date,
            )
            ClimaTimeSteppers.solve!(
                calendar_integrator(;
                    output_writer = writer,
                    output_short_name = "monthly_spin",
                    spinup = Dates.DateTime(2024, 2, 1),
                    pre_output_hook! = ClimaDiagnostics.average_pre_output_hook!,
                ),
            )
            close(writer)
            if ClimaComms.iamroot(context)
                NCDatasets.NCDataset(
                    joinpath(output_dir, "monthly_spin.nc"),
                ) do nc
                    # Only [Feb 1, Mar 1) is written, and its lower bound is Feb 1 (the window
                    # start), not the simulation start the old reconstruction would have used.
                    @test length(nc["time"]) == 1
                    @test nc["date_bnds"][:, 1] == [
                        Dates.DateTime(2024, 2, 1),
                        Dates.DateTime(2024, 3, 1),
                    ]
                    @test nc["date"][1] == Dates.DateTime(2024, 2, 1)
                    @test nc["time_bnds"][:, 1] == [31 * day, 60 * day]
                    @test nc["time"][1] == 31 * day
                end
            end
        end
    end

    @testset "instantaneous diagnostic keeps reconstructed bounds" begin
        # Instantaneous diagnostics must not get the accumulation window as bounds; they keep
        # the [previous, current] reconstruction.
        mktempdir() do output_dir
            output_dir = ClimaComms.bcast(context, output_dir)
            writer = ClimaDiagnostics.Writers.NetCDFWriter(
                ColumnCenterFiniteDifferenceSpace(),
                output_dir;
                num_points = (10,),
                start_date = const_start_date,
            )
            ClimaTimeSteppers.solve!(
                calendar_integrator(;
                    output_writer = writer,
                    output_short_name = "inst_spin",
                    spinup = Dates.DateTime(2024, 2, 1),
                    reduction_time_func = nothing,
                ),
            )
            close(writer)
            if ClimaComms.iamroot(context)
                NCDatasets.NCDataset(joinpath(output_dir, "inst_spin.nc")) do nc
                    # The t0 snapshot (written at construction) and the post-spinup Mar 1 one.
                    @test nc["time"][:] == [0.0, 60 * day]
                    @test nc["time_bnds"][:, 2] == [0.0, 60 * day]
                    @test nc["date_bnds"][:, 2] == [
                        Dates.DateTime(2024, 1, 1),
                        Dates.DateTime(2024, 3, 1),
                    ]
                end
            end
        end
    end
end
