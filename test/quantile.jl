using Test
import ClimaComms
import ClimaCore: Fields, Geometry, Spaces
import ClimaDiagnostics
import ClimaDiagnostics: ExactQuantile
import ClimaDiagnostics.Schedules
import ClimaDiagnostics.Writers
import ClimaUtilities.TimeManager: ITime
import Dates
import Random
import Statistics

include("TestTools.jl")

const context = ClimaComms.context()
const on_gpu = ClimaComms.device(context) isa ClimaComms.CUDADevice

# Statistics.quantile, with the choices ExactQuantile makes where it has none: NaN for a
# window with NaN or no samples, and the sample with all the weight when the other is
# infinite (Statistics computes 0 * Inf = NaN there)
function reference_quantile(x, q)
    (isempty(x) || any(isnan, x)) && return NaN
    v = sort(x)
    n = length(v)
    aleph = fma(n, q, 1 - q)
    j = clamp(trunc(Int, aleph), 1, n - 1)
    γ = clamp(aleph - j, 0.0, 1.0)
    a, b = n == 1 ? (v[1], v[1]) : (v[j], v[j + 1])
    if isfinite(a) && isfinite(b)
        return a ≈ b ? a + γ * (b - a) : (1 - γ) * a + γ * b
    end
    γ == 0 && return a + zero(γ)
    γ == 1 && return b + zero(γ)
    return (1 - γ) * a + γ * b
end

# The GPU fuses multiply-adds, so it can differ from the CPU in the last bit
matches(result, expected, scale) =
    isequal(result, expected) ||
    on_gpu && isapprox(result, expected; atol = 4eps(typeof(result)) * scale)

# Run the hooks over one window of samples, reusing `buffer` across windows
function window_quantile!(out, buffer, reduction, field, samples; name = "test")
    ClimaDiagnostics.reset_accumulator!(buffer, reduction)
    for (t, sample) in enumerate(samples)
        copyto!(parent(field), sample)
        ClimaDiagnostics.accumulate_sample!(buffer, reduction, field, t)
    end
    ClimaDiagnostics.finalize_accumulator!(
        out,
        buffer,
        reduction,
        length(samples),
        name,
    )
    return Array(parent(out))
end

# Quantile of the samples `xs` on a single point
function point_quantile(xs, q; batch_size = 0, max_samples = length(xs))
    FT = eltype(xs)
    space = Spaces.PointSpace(context, Geometry.ZPoint(zero(FT)))
    field = Fields.zeros(FT, space)
    reduction = ExactQuantile(q; batch_size)
    buffer = ClimaDiagnostics.allocate_accumulator(
        reduction,
        field,
        max(max_samples, 1),
        "test",
    )
    samples = [fill(x, size(parent(field))) for x in xs]
    return only(
        window_quantile!(similar(field), buffer, reduction, field, samples),
    )
end

const QS = (0.0, 0.05, 0.1, 0.5, 0.95, 0.999, 1.0)

@testset "Names and constructor" begin
    @test string(ExactQuantile(0.95)) == "p95"
    @test string(ExactQuantile(0.999)) == "p99p9"
    @test string(ExactQuantile(0.05)) == "p5"
    @test string(ExactQuantile(0.5)) == "p50"
    @test string(ExactQuantile(0)) == "p0"
    @test string(ExactQuantile(1)) == "p100"
    @test string(ExactQuantile(1 // 4)) == "p25"
    @test string(ExactQuantile(1e-9)) == "p0p0000001"
    @test string(ExactQuantile(1 - 1e-9)) == "p99p9999999"
    @test string(ExactQuantile(-0.0)) == "p0"
    @test string(ExactQuantile(0.95f0)) == "p95"
    @test ExactQuantile(0.95f0).q === 0.95

    @test ExactQuantile(0.95; batch_size = 16, max_samples = 3) ==
          ExactQuantile(0.95)
    @test ExactQuantile(0.95) != ExactQuantile(0.99)
    @test hash(ExactQuantile(0.95; batch_size = 16)) ==
          hash(ExactQuantile(0.95))

    @test_throws ArgumentError ExactQuantile(1.5)
    @test_throws ArgumentError ExactQuantile(NaN)
    @test_throws ArgumentError ExactQuantile(0.5; batch_size = -1)
    @test_throws ArgumentError ExactQuantile(0.5; max_samples = 0)
end

@testset "Reference matches Statistics" begin
    rng = Random.Xoshiro(1)
    for n in 1:40, q in QS
        x = randn(rng, n)
        @test reference_quantile(x, q) ≈ Statistics.quantile(x, q) atol =
            8eps() * maximum(abs, x)
    end
end

@testset "Hooks against the definition" begin
    rng = Random.Xoshiro(2)
    max_samples = 20
    for FT in (Float32, Float64)
        spaces = (
            Spaces.PointSpace(context, Geometry.ZPoint(zero(FT))),
            ColumnCenterFiniteDifferenceSpace(4, context; FT),
            Spaces.horizontal_space(
                SphericalShellSpace(;
                    nelements = 2,
                    zelem = 2,
                    npolynomial = 2,
                    context,
                    FT,
                ),
            ),
            SphericalShellSpace(;
                nelements = 2,
                zelem = 3,
                npolynomial = 2,
                context,
                FT,
            ),
        )
        for space in spaces, q in QS, batch_size in (0, 7)
            field = Fields.zeros(FT, space)
            out = similar(field)
            reduction = ExactQuantile(q; batch_size)
            buffer = ClimaDiagnostics.allocate_accumulator(
                reduction,
                field,
                max_samples,
                "test",
            )
            # Windows shorter than the heap, longer than it, and at the limit
            for n in (1, 2, 9, max_samples)
                samples = [randn(rng, FT, size(parent(field))) for _ in 1:n]
                result =
                    window_quantile!(out, buffer, reduction, field, samples)
                expected = map(eachindex(result)) do i
                    FT(reference_quantile([s[i] for s in samples], q))
                end
                scale = maximum(s -> maximum(abs, s), samples)
                @test all(matches.(vec(result), expected, scale))
            end
        end
    end
end

@testset "Edge cases" begin
    cases(FT, big = floatmax(FT)) = (
        # The sample with all the weight wins over an infinite neighbour
        (FT[1, 2, Inf], 0.5, 2),
        (FT[1, Inf], 0.0, 1),
        (FT[-Inf, 1], 1.0, 1),
        (FT[Inf, Inf], 0.0, Inf),
        (FT[1, Inf], 1.0, Inf),
        (FT[1, 2, Inf], 0.75, Inf),
        (FT[-Inf, -Inf, 0.1], 0.1, -Inf),
        # No overflow between huge values of opposite sign
        (FT[-big, big], 0.5, 0),
        (FT[-big, -big, big], 0.5, -big),
        # Zeros stay positive below the median
        (zeros(FT, 10), 0.0, 0),
        (zeros(FT, 10), 0.1, 0),
        (zeros(FT, 10), 0.9, 0),
        # NaN anywhere gives NaN
        (FT[NaN, 1, 2], 0.5, NaN),
        (FT[1, NaN, 2], 0.1, NaN),
        (FT[1, 2, 3, NaN], 0.95, NaN),
        # One sample, and ties
        (FT[3], 0.3, 3),
        (FT[1, 1, 2, 2, 2, 3], 0.5, 2),
        (FT[1, 1, 2, 2, 2, 3], 0.1, 1),
    )
    for FT in (Float32, Float64), (xs, q, expected) in cases(FT)
        for batch_size in (0, 2)
            @test isequal(point_quantile(xs, q; batch_size), FT(expected))
        end
    end
    # A window with no samples
    @test isnan(point_quantile(Float64[], 0.5; max_samples = 1))
    # p0 and p100 are the min and max
    rng = Random.Xoshiro(3)
    xs = randn(rng, 50)
    @test point_quantile(xs, 0.0) == minimum(xs)
    @test point_quantile(xs, 1.0) == maximum(xs)
end

@testset "Windows longer than the buffer" begin
    rng = Random.Xoshiro(4)
    space = Spaces.PointSpace(context, Geometry.ZPoint(0.0))
    field = Fields.zeros(space)
    for q in (0.1, 0.5)
        reduction = ExactQuantile(q)
        buffer = ClimaDiagnostics.allocate_accumulator(reduction, field, 10, "")
        n = buffer.capacity
        for num_samples in 11:16
            xs = randn(rng, num_samples)
            samples = [fill(x, size(parent(field))) for x in xs]
            needed =
                ClimaDiagnostics._num_slots_needed(q, num_samples, buffer.keep)
            shortfall = needed - n
            name = "x_p$(q)_$num_samples"
            run() = only(
                window_quantile!(
                    similar(field),
                    buffer,
                    reduction,
                    field,
                    samples;
                    name,
                ),
            )
            if shortfall <= 0
                @test run() == reference_quantile(xs, q)
            elseif shortfall <= 2
                # The closest value the buffer holds is its n-th most extreme sample
                closest = sort(xs; rev = q >= 0.5)[n]
                @test (@test_logs (
                    :warn,
                    Regex("^$name: .*short by $shortfall"),
                ) run()) == closest
            else
                @test isnan(@test_logs (:warn, Regex("^$name: .*NaN")) run())
            end
        end
    end
end

@testset "Pointwise heaps keep the most extreme samples" begin
    rng = Random.Xoshiro(7)
    num_points = 3
    capacity = 8
    keeps = (ClimaDiagnostics.Largest(), ClimaDiagnostics.Smallest())
    for keep in keeps, batch_capacity in (0, 3), n in (4, 8, 20)
        heaps = ClimaDiagnostics.PointwiseHeaps(
            keep,
            zeros(num_points),
            capacity;
            batch_capacity,
            device = ClimaComms.CPUSingleThreaded(),
        )
        samples = [randn(rng, num_points) for _ in 1:n]
        foreach(sample -> push!(heaps, sample), samples)
        @test length(heaps) == min(n, capacity)
        for m in (min(n, capacity), 3, 1)
            ClimaDiagnostics.pop_until!(heaps, m)
            roots, nexts = zeros(num_points), zeros(num_points)
            ClimaDiagnostics.map_top_two!((root, _) -> root, roots, heaps)
            ClimaDiagnostics.map_top_two!((_, next) -> next, nexts, heaps)
            for i in 1:num_points
                most_extreme_first = sort(
                    [s[i] for s in samples];
                    rev = keep isa ClimaDiagnostics.Largest,
                )
                @test roots[i] == most_extreme_first[m]
                @test nexts[i] == most_extreme_first[max(m - 1, 1)]
            end
        end
        # Emptying starts a new window
        empty!(heaps)
        @test length(heaps) == 0
    end
end

@testset "Plain arrays" begin
    storage = [0.0]
    reduction = ExactQuantile(0.5)
    buffer =
        ClimaDiagnostics.allocate_accumulator(reduction, storage, 5, "plain")
    xs = [3.0, 1.0, 4.0, 1.0, 5.0]
    for (t, x) in enumerate(xs)
        storage .= x
        ClimaDiagnostics.accumulate_sample!(buffer, reduction, storage, t)
    end
    ClimaDiagnostics.finalize_accumulator!(storage, buffer, reduction, 5, "")
    @test storage[] == 3.0
end

@testset "Samples per window from the schedules" begin
    every_month = Schedules.EveryCalendarDtSchedule(
        Dates.Month(1);
        start_date = Dates.DateTime(2000),
    )
    every_step = Schedules.EveryStepSchedule()
    bound = Schedules.max_samples_per_window
    @test bound(
        Schedules.DivisorSchedule(2),
        Schedules.DivisorSchedule(10),
        nothing,
    ) == 5
    @test bound(every_step, Schedules.DivisorSchedule(10), nothing) == 10
    @test bound(
        Schedules.EveryDtSchedule(3600.0),
        Schedules.DivisorSchedule(10),
        nothing,
    ) == 10
    @test bound(Schedules.EveryDtSchedule(3600.0), every_month, nothing) ==
          744 + 2
    @test bound(every_step, every_month, 450.0) == 5952 + 2
    @test isnothing(bound(every_step, every_month, nothing))
    @test bound(
        Schedules.EveryDtSchedule(ITime(3600)),
        Schedules.EveryDtSchedule(ITime(86400)),
        ITime(450),
    ) == 24 + 2
    @test isnothing(bound(integrator -> true, every_month, 450.0))
end

@testset "Through the handler" begin
    dt = 1.0
    num_steps = 30
    space = ColumnCenterFiniteDifferenceSpace(4, context)
    rng = Random.Xoshiro(5)
    values =
        [randn(rng, size(parent(Fields.zeros(space)))) for _ in 0:num_steps]
    value_at(t) = values[round(Int, t / dt) + 1]
    function compute!(out, u, p, t)
        if isnothing(out)
            out = Fields.zeros(space)
            copyto!(parent(out), value_at(t))
            return out
        end
        copyto!(parent(out), value_at(t))
        return nothing
    end
    variable = ClimaDiagnostics.DiagnosticVariable(; compute!, short_name = "x")
    # A diagnostic that returns a plain array
    plain_variable = ClimaDiagnostics.DiagnosticVariable(;
        compute! = (out, u, p, t) ->
            isnothing(out) ? [value_at(t)[1]] :
            (out .= value_at(t)[1]; nothing),
        short_name = "y",
    )

    writer = Writers.DictWriter()
    diagnostic(variable, q, compute, output; batch_size = 0) =
        ClimaDiagnostics.ScheduledDiagnostic(;
            variable,
            output_writer = writer,
            reduction_time_func = ExactQuantile(q; batch_size),
            compute_schedule_func = compute,
            output_schedule_func = output,
        )
    diagnostics = [
        diagnostic(
            variable,
            0.95,
            Schedules.EveryStepSchedule(),
            Schedules.DivisorSchedule(10),
        ),
        diagnostic(
            variable,
            0.1,
            Schedules.DivisorSchedule(2),
            Schedules.DivisorSchedule(10);
            batch_size = 3,
        ),
        # Step-based compute with time-based output, sized from dt
        diagnostic(
            variable,
            0.5,
            Schedules.EveryStepSchedule(),
            Schedules.EveryDtSchedule(7.0),
        ),
        diagnostic(
            plain_variable,
            0.5,
            Schedules.EveryStepSchedule(),
            Schedules.DivisorSchedule(5),
        ),
        # Same q as the first one, so a duplicate
        diagnostic(
            variable,
            0.95,
            Schedules.EveryStepSchedule(),
            Schedules.DivisorSchedule(10);
            batch_size = 16,
        ),
    ]
    Y = Fields.FieldVector(; my_var = Fields.zeros(space))
    handler = @test_logs (
        :warn,
        "Given list of diagnostics contains duplicates, removing them",
    ) ClimaDiagnostics.DiagnosticsHandler(diagnostics, Y, (;), 0.0; dt)
    @test length(handler.scheduled_diagnostics) == 4

    for step in 1:num_steps
        integrator = (; u = Y, p = (;), t = step * dt, step, dt)
        ClimaDiagnostics.orchestrate_diagnostics(integrator, handler)
    end

    # Each window covers the steps in (start, end]
    function check(diag, steps_per_window, compute_every)
        results = writer[diag.output_short_name]
        @test length(results) == num_steps ÷ steps_per_window
        for (start, result) in results
            first_step = round(Int, start / dt)
            steps =
                (first_step + compute_every):compute_every:(first_step + steps_per_window)
            q = diag.reduction_time_func.q
            got = Array(result isa Fields.Field ? parent(result) : result)
            for i in eachindex(got)
                expected =
                    reference_quantile([values[s + 1][i] for s in steps], q)
                @test matches(got[i], expected, 10.0)
            end
        end
    end
    p95, p10, p50, plain_p50 = handler.scheduled_diagnostics
    check(p95, 10, 1)
    check(p10, 10, 2)
    check(p50, 7, 1)
    check(plain_p50, 5, 1)
    # The buffers are sized from the schedules
    needed =
        ClimaDiagnostics._num_slots_needed(0.95, 10, ClimaDiagnostics.Largest())
    @test handler.accumulators[1].capacity == needed
end

@testset "Errors at construction" begin
    space = ColumnCenterFiniteDifferenceSpace(4, context)
    variable = ClimaDiagnostics.DiagnosticVariable(;
        compute = (u, p, t) -> Fields.zeros(space),
        short_name = "x",
    )
    Y = Fields.FieldVector(; my_var = Fields.zeros(space))
    diagnostic(compute, q) = ClimaDiagnostics.ScheduledDiagnostic(;
        variable,
        output_writer = Writers.DictWriter(),
        reduction_time_func = q,
        compute_schedule_func = compute,
        output_schedule_func = Schedules.EveryDtSchedule(10.0),
    )
    # The schedules don't bound the samples per window
    @test_throws ErrorException ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic(integrator -> true, ExactQuantile(0.95))],
        Y,
        (;),
        0.0;
        dt = 1.0,
    )
    @test_throws ErrorException ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic(Schedules.EveryStepSchedule(), ExactQuantile(0.95))],
        Y,
        (;),
        0.0,
    )
    # Unless max_samples is given
    handler = ClimaDiagnostics.DiagnosticsHandler(
        [diagnostic(integrator -> true, ExactQuantile(0.95; max_samples = 10))],
        Y,
        (;),
        0.0;
        dt = 1.0,
    )
    @test length(handler.scheduled_diagnostics) == 1
end
