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

@testset "Compute fields" begin
    var_compute_in_place = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        standard_name = "my_test",
        units = "m",
        comments = "It works!",
        compute! = (out, u, p, t) -> isnothing(out) ? 1 : (out .= 1),
    )
    var_compute = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        standard_name = "my_test",
        units = "m",
        comments = "It works!",
        compute = (u, p, t) -> 1,
    )
    var_lazy_compute = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        standard_name = "my_test",
        units = "m",
        comments = "It works!",
        compute = (u, p, t) -> Base.Broadcast.broadcasted(+, 0, 1),
    )

    function create_diag(variable)
        return ScheduledDiagnostics.ScheduledDiagnostic(;
            variable = variable,
            output_writer = Writers.DummyWriter(),
            output_schedule_func = Schedules.EveryDtSchedule(5),
            compute_schedule_func = Schedules.EveryDtSchedule(5),
        )
    end

    u = nothing
    p = nothing
    t = nothing
    vars = (var_compute_in_place, var_compute, var_lazy_compute)

    for var in vars
        # Test allocating version
        @test isone(ClimaDiagnostics.compute_field(create_diag(var), u, p, t))

        # Test non allocating version
        dest = [0]
        ClimaDiagnostics.compute_field!(dest, create_diag(var), u, p, t)
        @test isone(first(dest))
    end
end
