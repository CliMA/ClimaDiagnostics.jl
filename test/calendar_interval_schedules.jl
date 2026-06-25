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

    # Step day-by-day as the orchestration does (`should_accumulate` then the schedule call),
    # collecting the closed windows and the per-step accumulation flag.
    function drive(schedule; ndays)
        fired = Tuple{Float64, Any, Bool}[]
        accumulate_log = Bool[]
        for k in 0:ndays
            integrator = (; t = k * day)
            acc = should_accumulate(schedule, integrator)
            push!(accumulate_log, acc)
            schedule(integrator) &&
                push!(fired, (k * day, output_period_bounds(schedule), acc))
        end
        return fired, accumulate_log
    end

    @testset "monthly windows (no spinup)" begin
        schedule = CalendarIntervalSchedule(Dates.Month(1); start_date)
        fired, acc = drive(schedule; ndays = 70)

        @test length(fired) == 2
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
        @test all(acc)                                       # every window active
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
        # spinup = Feb 1 makes (Jan 1, Feb 1] inactive; the first output window is (Feb 1, Mar 1].
        schedule = CalendarIntervalSchedule(
            Dates.Month(1);
            start_date,
            spinup = Dates.DateTime(2024, 2, 1),
        )
        fired, acc = drive(schedule; ndays = 70)

        @test length(fired) == 1
        @test fired[1][1] == 60 * day
        @test fired[1][2] == ScheduleInterval(
            Dates.DateTime(2024, 2, 1),
            Dates.DateTime(2024, 3, 1),
        )
        # `should_accumulate` runs before the schedule advances, so day 31 (the boundary sample
        # of the inactive January window) is also excluded; accumulation starts on day 32.
        @test all(.!acc[1:32])
        @test all(acc[33:end])
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
end

# End-to-end: a sample is folded only when the compute schedule fires AND the time is in an
# active window. With `value = t -> 1` and a `+` reduction (no hook), the output equals the
# number of folded samples. dt = 1 day, so boundaries land on steps 31 (Feb 1) and 60 (Mar 1).
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
            output_schedule_func = CalendarIntervalSchedule(
                Dates.Month(1);
                start_date = const_start_date,
                spinup,
            ),
            pre_output_hook!,
            output_short_name,
        )
        return ClimaDiagnostics.IntegratorWithDiagnostics(
            ClimaTimeSteppers.init(args...; kwargs...),
            [diag],
        )
    end

    @testset "sample gating (sum reduction)" begin
        # No spinup: (Jan 1, Feb 1] holds the t0 seed + 31 daily samples (32); (Feb 1, Mar 1] holds 29.
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
            @test parent(out[31 * day])[1] == 32
            @test parent(out[60 * day])[1] == 29
        end

        # spinup = Feb 1: January never accumulates, so only (Feb 1, Mar 1] is written (29 samples).
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

    @testset "non-additive reduction with spinup (min)" begin
        # Each sample equals `t`; `min` over (Feb 1, Mar 1] is day 32. A bad identity seed or a
        # leaked pre-spinup sample (e.g. the t0 = 0 seed) would drag the minimum to 0.
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
            @test parent(out[60 * day])[1] == 32 * day
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
                    # Only (Feb 1, Mar 1] is written, and its lower bound is Feb 1 (the window
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
