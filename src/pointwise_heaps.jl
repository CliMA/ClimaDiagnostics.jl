# One binary heap per point of an array, keeping the `capacity` largest or smallest samples
# of a stream of arrays. The heap of the point with linear index i is slots[i, 1:capacity],
# so each slot is a contiguous array and GPU memory access stays coalesced. The root
# slots[i, 1] is the least extreme sample kept, so a new sample enters only if it is more
# extreme than the root, which it replaces.
#
# Per window: `empty!`, `push!` each sample, then `pop_until!(heaps, m)` and `map_top_two!`,
# which can repeat with decreasing m.
#
# Samples are copied into the free slots after the heap, the batch. Once all
# `capacity + batch_capacity` slots are full, the batch is pushed onto the heaps, which are
# built first if needed. With `batch_capacity = 0`, later samples are pushed one at a time.
# `pop_until!` pushes the last, partial batch.
#
# The functions without a leading underscore act on all points at once. Those ending in
# `_point` act on the heap of one point, inside a kernel.

# Which samples the heaps keep
struct Largest end
struct Smallest end

struct PointwiseHeaps{K <: Union{Largest, Smallest}, A <: AbstractArray, D}
    keep::K
    slots::A
    capacity::Int
    batch_capacity::Int
    device::D
    # slots[i, 1:heap_size] form the heap of point i, and the next batch_size slots hold
    # samples copied in but not yet pushed onto it
    heap_size::Base.RefValue{Int}
    batch_size::Base.RefValue{Int}
end

function PointwiseHeaps(keep, data, capacity; batch_capacity = 0, device)
    slots = similar(data, length(data), capacity + batch_capacity)
    return PointwiseHeaps(
        keep,
        slots,
        capacity,
        batch_capacity,
        device,
        Ref(0),
        Ref(0),
    )
end

function Base.empty!(heaps::PointwiseHeaps)
    heaps.heap_size[] = 0
    heaps.batch_size[] = 0
    return heaps
end

# Number of samples each heap holds once the batch is pushed
Base.length(heaps::PointwiseHeaps) =
    min(heaps.heap_size[] + heaps.batch_size[], heaps.capacity)

# Copying samples into free slots and pushing them once the slots are full keeps the sample
# count out of the kernels, since CUDA caches a kernel for each distinct value they capture
function Base.push!(heaps::PointwiseHeaps, data)
    (; slots, capacity, batch_capacity) = heaps
    h, b = heaps.heap_size[], heaps.batch_size[]
    if h + b < capacity + batch_capacity
        _copy_to_slot!(slots, data, h + b + 1)
        heaps.batch_size[] = b + 1
        h + b + 1 == capacity + batch_capacity && pop_until!(heaps, capacity)
    else
        _pushpop!(heaps, data)
    end
    return heaps
end

# Push the batch and pop roots until each heap holds its m most extreme samples, where m is
# at most `length(heaps)`. Each root is then the m-th most extreme sample of its point.
function pop_until!(heaps::PointwiseHeaps, m)
    (; keep, slots, device) = heaps
    heaps.heap_size[] == 0 && _heapify!(heaps, m)
    h, b = heaps.heap_size[], heaps.batch_size[]
    h == m && b == 0 && return heaps
    ClimaComms.@threaded device for i in _points(heaps)
        _pop_until_point!(slots, i, h, b, m, keep)
    end
    heaps.heap_size[] = m
    heaps.batch_size[] = 0
    return heaps
end

# After `pop_until!`, set out[i] = f(root, next root) for the heap of each point i
function map_top_two!(f, out, heaps::PointwiseHeaps)
    (; keep, slots, device) = heaps
    m = length(heaps)
    ClimaComms.@threaded device for i in _points(heaps)
        root, next = _top_two_point(slots, i, m, keep)
        @inbounds out[i] = f(root, next)
    end
    return out
end

# The heap order: whether x is less extreme than y, so belongs nearer the root
@inline _compare(::Largest, x, y) = x < y
@inline _compare(::Smallest, x, y) = x > y

# A linear index saves each GPU thread converting it to a Cartesian index and back
_points(heaps) = axes(heaps.slots, 1)

_copy_to_slot!(slots, data, slot) =
    copyto!(slots, (slot - 1) * length(data) + 1, data, 1, length(data))

# Build heaps from the first m samples of the batch
function _heapify!(heaps, m)
    (; keep, slots, device) = heaps
    ClimaComms.@threaded device for i in _points(heaps)
        _heapify_point!(slots, i, m, keep)
    end
    heaps.heap_size[] = m
    heaps.batch_size[] -= m
    return nothing
end

function _pushpop!(heaps, data)
    (; keep, slots, capacity, device) = heaps
    ClimaComms.@threaded device for i in _points(heaps)
        @inbounds _pushpop_point!(slots, i, data[i], capacity, keep)
    end
    return nothing
end

# Push `value` onto the full heap slots[i, 1:n] and pop the least extreme sample, as Python's
# heappushpop. A NaN turns the whole heap to NaN, so any result read from it is NaN.
@inline function _pushpop_point!(slots, i, value, n, keep)
    @inbounds if isnan(value)
        if !isnan(slots[i, 1])
            for k in 1:n
                slots[i, k] = value
            end
        end
    elseif _compare(keep, slots[i, 1], value)
        _sift_down_point!(slots, i, 1, value, n, keep)
    end
    return nothing
end

@inline function _heapify_point!(slots, i, n, keep)
    has_nan = false
    @inbounds for k in 1:n
        has_nan |= isnan(slots[i, k])
    end
    if has_nan
        @inbounds for k in 1:n
            slots[i, k] = NaN
        end
    else
        @inbounds for k in (n ÷ 2):-1:1
            _sift_down_point!(slots, i, k, slots[i, k], n, keep)
        end
    end
    return nothing
end

# Put `value` at slot k of the heap slots[i, 1:n] and move it down to its place
@inline function _sift_down_point!(slots, i, k, value, n, keep)
    @inbounds while true
        child = 2k
        child > n && break
        if child < n && _compare(keep, slots[i, child + 1], slots[i, child])
            child += 1
        end
        _compare(keep, slots[i, child], value) || break
        slots[i, k] = slots[i, child]
        k = child
    end
    @inbounds slots[i, k] = value
    return nothing
end

# Pop the heap slots[i, 1:h] down to its m most extreme samples, then push the batch
# slots[i, h .+ (1:b)] onto it
@inline function _pop_until_point!(slots, i, h, b, m, keep)
    @inbounds for k in h:-1:(m + 1)
        _sift_down_point!(slots, i, 1, slots[i, k], k - 1, keep)
    end
    @inbounds for s in (h + 1):(h + b)
        _pushpop_point!(slots, i, slots[i, s], m, keep)
    end
    return nothing
end

# The root of the heap slots[i, 1:m] and the next root
@inline function _top_two_point(slots, i, m, keep)
    @inbounds root = slots[i, 1]
    m == 1 && return root, root
    @inbounds next = slots[i, 2]
    if m > 2 && @inbounds _compare(keep, slots[i, 3], next)
        @inbounds next = slots[i, 3]
    end
    return root, next
end
