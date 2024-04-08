using Test
import SciMLBase

import ClimaDiagnostics.Callbacks:
    Callback, CallbackOrchestrator, DivisorSchedule, EveryDtSchedule

include("TestTools.jl")

@testset "Callback" begin
    # Test a callback that is called at every iteration by always returning true

    called = Ref(0)
    function callback_func(integrator)
        called[] += 1
    end

    scheduled_func = (_integrator) -> true
    callback_everystep = Callback(callback_func, scheduled_func)

    t0 = 0.0
    tf = 1.0
    dt = 1e-3

    space = ColumnCenterFiniteDifferenceSpace()
    args, kwargs = create_problem(space; t0, tf, dt)

    expected_called = convert(Int, (tf - t0) / dt)

    cb = CallbackOrchestrator([callback_everystep])
    SciMLBase.solve(args...; kwargs..., callback = cb)

    @test called[] == expected_called
end

@testset "Schedules" begin
    # Test a callback that is called at every other iteration with DivisorSchedule

    called = Ref(0)
    function callback_func0(integrator)
        called[] += 1
    end

    divisor = 2
    scheduled_func = DivisorSchedule(divisor)
    @test "$scheduled_func" == "2it"

    callback_everystep = Callback(callback_func0, scheduled_func)

    t0 = 0.0
    tf = 1.0
    dt = 1e-3

    space = ColumnCenterFiniteDifferenceSpace()
    args, kwargs = create_problem(space; t0, tf, dt)

    expected_called = convert(Int, (tf - t0) / (divisor * dt))

    cb = CallbackOrchestrator([callback_everystep])
    SciMLBase.solve(args...; kwargs..., callback = cb)

    @test called[] == expected_called

    # EveryDtSchedule

    called = Ref(0)
    function callback_func(integrator)
        called[] += 1
    end

    called2 = Ref(0)
    function callback_func2(integrator)
        called2[] += 1
    end

    dt_callback = 0.2
    scheduled_func = EveryDtSchedule(dt_callback)
    @test "$scheduled_func" == "0.2s"

    dt_callback2 = 0.3
    t_start2 = 0.1
    scheduled_func2 = EveryDtSchedule(dt_callback2; t_start = t_start2)

    callback_dt = Callback(callback_func, scheduled_func)
    callback_dt2 = Callback(callback_func2, scheduled_func2)

    args, kwargs = create_problem(space; t0, tf, dt)

    expected_called = convert(Int, (tf - t0) / dt_callback)
    expected_called2 = convert(Int, floor((tf - t0 - t_start2) / dt_callback2))

    cb = CallbackOrchestrator([callback_dt, callback_dt2])

    SciMLBase.solve(args...; kwargs..., callback = cb)

    @test called[] == expected_called
    @test called2[] == expected_called2
end
