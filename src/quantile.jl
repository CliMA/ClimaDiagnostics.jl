# Exact q-quantile over each output window, as defined by `Statistics.quantile`. Each point
# keeps its most extreme samples in `PointwiseHeaps`, sized for the most samples a window
# can have: `max_samples` if given, else the bound the handler derives from the schedules.
# A positive `batch_size` inserts samples in batches, which is faster on GPUs but costs
# `batch_size` more copies of the field.
struct ExactQuantile
    q::Float64
    max_samples::Union{Int, Nothing}
    batch_size::Int
    function ExactQuantile(q::Real; max_samples = nothing, batch_size = 0)
        q = _to_float64(q)
        0 <= q <= 1 || throw(ArgumentError("q must be in [0, 1], got $q"))
        something(max_samples, 1) >= 1 || throw(
            ArgumentError("max_samples must be positive, got $max_samples"),
        )
        batch_size >= 0 || throw(
            ArgumentError("batch_size must be non-negative, got $batch_size"),
        )
        return new(q + 0.0, max_samples, batch_size)
    end
end

# A float is read as the decimal it prints as, so that 0.95f0 means 0.95
_to_float64(q::AbstractFloat) = parse(Float64, string(q))
_to_float64(q::Real) = Float64(q)

# The other fields only change how the result is computed, not the result
Base.:(==)(r1::ExactQuantile, r2::ExactQuantile) = r1.q == r2.q
Base.hash(r::ExactQuantile, h::UInt) = hash(r.q, hash(:ExactQuantile, h))

# The name is 100q with "p" for the decimal point, e.g. p95, p99p9 and p0p0000001 (q = 1e-9).
# It is built from the shortest decimal that identifies q, so distinct q get distinct names.
function Base.show(io::IO, reduction::ExactQuantile)
    str = string(reduction.q)
    mantissa, power = occursin('e', str) ? split(str, 'e') : (str, "0")
    whole, fraction = split(mantissa, '.')
    all_digits = whole * fraction
    point = length(whole) + parse(Int, power) + 2
    if point < 1
        all_digits = "0"^(1 - point) * all_digits
        point = 1
    end
    all_digits = rpad(all_digits, point, '0')
    integer_part = lstrip(all_digits[1:point], '0')
    fraction_part = rstrip(all_digits[(point + 1):end], '0')
    print(io, "p", isempty(integer_part) ? "0" : integer_part)
    isempty(fraction_part) || print(io, "p", fraction_part)
end

# Below the median, the heaps keep the smallest samples, otherwise the largest
_samples_to_keep(reduction::ExactQuantile) =
    reduction.q < 0.5 ? Smallest() : Largest()

# The q-quantile of n samples interpolates with weight γ between the j-th and (j + 1)-th
# smallest samples, as in Statistics.quantile
function _interpolation_point(q, n)
    aleph = fma(n, q, 1 - q)
    j = clamp(trunc(Int, aleph), 1, n - 1)
    return j, clamp(aleph - j, 0.0, 1.0)
end

# Number of most extreme samples a heap needs for the two samples around the quantile to be
# its least extreme ones
function _num_slots_needed(q, n, keep)
    n == 1 && return 1
    j, _ = _interpolation_point(q, n)
    return keep isa Smallest ? j + 1 : n - j + 1
end

# Plain-array diagnostics have no ClimaComms device of their own
_quantile_device(storage) = ClimaComms.device(storage)
_quantile_device(::Array) = ClimaComms.CPUSingleThreaded()

function allocate_accumulator(
    reduction::ExactQuantile,
    storage,
    max_samples,
    name,
)
    max_samples = something(reduction.max_samples, max_samples, Some(nothing))
    isnothing(max_samples) && error(
        "$name: cannot tell from its schedules how many samples an output window can \
         have. Pass dt to the DiagnosticsHandler, or max_samples to ExactQuantile",
    )
    keep = _samples_to_keep(reduction)
    capacity =
        maximum(n -> _num_slots_needed(reduction.q, n, keep), 1:max_samples)
    return PointwiseHeaps(
        keep,
        parent(storage),
        capacity;
        batch_capacity = reduction.batch_size,
        device = _quantile_device(storage),
    )
end

function reset_accumulator!(heaps::PointwiseHeaps, ::ExactQuantile)
    empty!(heaps)
    return nothing
end

function accumulate_sample!(heaps::PointwiseHeaps, ::ExactQuantile, storage, _)
    push!(heaps, parent(storage))
    return nothing
end

# Statistics.quantile's interpolation, except that the sample with all the weight is
# returned when the other one is infinite, where Statistics computes 0 * Inf = NaN
@inline function _interpolate(a, b, γ)
    if isfinite(a) && isfinite(b)
        return a ≈ b ? a + γ * (b - a) : (1 - γ) * a + γ * b
    end
    γ == 0 && return a + zero(γ)
    γ == 1 && return b + zero(γ)
    return (1 - γ) * a + γ * b
end

# A window whose heaps are short by at most this many samples gets the closest value they
# hold; one short by more gets NaN
const MAX_SHORTFALL = 2

function finalize_accumulator!(
    storage,
    heaps::PointwiseHeaps,
    reduction::ExactQuantile,
    num_samples,
    name,
)
    out = parent(storage)
    if num_samples == 0
        fill!(out, NaN)
        return nothing
    end
    available = length(heaps)
    needed = _num_slots_needed(reduction.q, num_samples, heaps.keep)
    shortfall = needed - available
    if shortfall > 0
        approximate = shortfall <= MAX_SHORTFALL
        result =
            approximate ? "the closest value it holds (short by $shortfall)" :
            "NaN"
        message = "$name: a window had $num_samples samples, more than \
                   $reduction can hold exactly, so it wrote $result. Check that \
                   the dt given to the DiagnosticsHandler matches the \
                   integrator, or pass a larger max_samples. This warning is \
                   shown once."
        @warn message maxlog = 1 _id = (name, approximate)
        if !approximate
            fill!(out, NaN)
            return nothing
        end
    end
    pop_until!(heaps, min(needed, available))
    if shortfall > 0
        map_top_two!((root, _) -> root, out, heaps)
    else
        _, γ = _interpolation_point(reduction.q, num_samples)
        map_top_two!(out, heaps) do root, next
            # The j-th and (j + 1)-th smallest samples
            a, b = minmax(root, next)
            return _interpolate(a, b, γ)
        end
    end
    return nothing
end
