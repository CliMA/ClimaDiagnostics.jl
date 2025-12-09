using Test

import ClimaDiagnostics
import ClimaDiagnostics:
    DiagnosticVariables, Writers, ScheduledDiagnostics, Schedules

@testset "Check dt schedules" begin
    dt = nothing
    @test isnothing(ClimaDiagnostics._check_dt_schedules(dt, nothing))

    var = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        standard_name = "my_test",
        units = "m",
        comments = "It works!",
        compute! = (out, u, p, t) -> 1,
    )

    function create_diag(output_schedule_func, compute_schedule_func)
        return ScheduledDiagnostics.ScheduledDiagnostic(;
            variable = var,
            output_writer = Writers.DummyWriter(),
            output_schedule_func = output_schedule_func,
            compute_schedule_func = compute_schedule_func,
        )
    end
    dt = 5
    output_schedule_func = Schedules.EveryDtSchedule(5)
    compute_schedule_func = Schedules.EveryDtSchedule(5)
    diag = create_diag(output_schedule_func, compute_schedule_func)
    @test isnothing(ClimaDiagnostics._check_dt_schedules(dt, [diag]))

    schedule1 = Schedules.EveryDtSchedule(5)
    schedule2 = Schedules.EveryDtSchedule(11)
    for (output_schedule_func, compute_schedule_func) in
        ((schedule1, schedule2), (schedule2, schedule1))
        diag = create_diag(output_schedule_func, compute_schedule_func)
        @test_throws ErrorException ClimaDiagnostics._check_dt_schedules(
            dt,
            [diag],
        )
    end
end
