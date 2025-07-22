using StaticArrays

# This uses fixed-size histograms to track value distributions

"""
    BinnedAccumulator{T, N_BINS, N_EDGES}

A memory-efficient, `isbits` accumulator that maintains a histogram of values.
The number of bins (`N_BINS`) and bin edges (`N_EDGES`) are type parameters
to allow for static allocation. `N_EDGES` must equal `N_BINS + 1`.
"""
struct BinnedAccumulator{T, N_BINS, N_EDGES}
    bin_edges::SVector{N_EDGES, T}
    bin_counts::SVector{N_BINS, T}
    total_count::T
end

# Make BinnedAccumulator broadcast like a scalar
Base.Broadcast.broadcastable(bda::BinnedAccumulator) = Ref(bda)

# Constructor for empty binned accumulator from a standard Vector
function BinnedAccumulator(bin_edges::Vector{T}) where {T}
    n_edges = length(bin_edges)
    n_bins = n_edges - 1
    return BinnedAccumulator{T, n_bins, n_edges}(
        SVector{n_edges}(bin_edges),
        zeros(SVector{n_bins, T}),
        zero(T),
    )
end

# Constructor for resetting the accumulator (used by pre_output_hook!)
function BinnedAccumulator(bin_edges::SVector{N_EDGES, T}) where {N_EDGES, T}
    n_bins = N_EDGES - 1
    return BinnedAccumulator{T, n_bins, N_EDGES}(
        bin_edges,
        zeros(SVector{n_bins, T}),
        zero(T),
    )
end

# Copy constructor
function Base.copy(bda::BinnedAccumulator{T, NB, NE}) where {T, NB, NE}
    return BinnedAccumulator{T, NB, NE}(
        bda.bin_edges,
        bda.bin_counts,
        bda.total_count,
    )
end

function Base.show(
    io::IO,
    bda::BinnedAccumulator{T, N_BINS, NE},
) where {T, N_BINS, NE}
    print(io, "BinnedAccumulator{$T, $N_BINS bins}($(bda.total_count) values)")
end

"""
    BinnedReducer

A callable struct that acts as a reduction function for binned diagnostics.
It holds the initial accumulator state, which is necessary for the diagnostics
system to initialize and reset the diagnostic correctly.
"""
struct BinnedReducer{Acc <: BinnedAccumulator}
    initial_accumulator::Acc
end

# This makes instances of `BinnedReducer` callable.
(reducer::BinnedReducer)(accumulator, new_value) =
    binned_reduction(accumulator, new_value)

# This is the key: we define how to get the identity *from an instance*
# of the reducer, which resolves the scoping issue.
identity_of_reduction(reducer::BinnedReducer) = reducer.initial_accumulator

"""
    binned_reduction(accumulator::BinnedAccumulator, new_value)

The reduction function that increments the appropriate histogram bin.
This function is called element-wise across the diagnostic field.
"""
function binned_reduction(
    accumulator::BinnedAccumulator{T, N_BINS, N_EDGES},
    new_value::T,
) where {T, N_BINS, N_EDGES}
    # Find which bin the new value belongs to
    bin_idx = bin_index(accumulator.bin_edges, new_value)
    isnothing(bin_idx) && return accumulator

    # SVectors are immutable, so we create a new one with the updated count.
    new_counts = setindex(
        accumulator.bin_counts,
        accumulator.bin_counts[bin_idx] + one(T),
        bin_idx,
    )

    return BinnedAccumulator{T, N_BINS, N_EDGES}(
        accumulator.bin_edges,
        new_counts,
        accumulator.total_count + one(T),
    )
end

"""
    bin_index(bin_edges, value)

Find which bin a value belongs to. Returns bin index (1-based).
Bins are left-inclusive: `[edge, next_edge)`.
Values outside the bin range are ignored.
"""
function bin_index(bin_edges::SVector{N, T}, value::T) where {N, T}
    if value < first(bin_edges) || value >= last(bin_edges)
        return nothing
    end
    # searchsortedlast returns the index of the last value in `bin_edges` that is <= `value`.
    return searchsortedlast(bin_edges, value)
end

"""
    binned_diagnostic(
        variable::DiagnosticVariable,
        bin_edges::Vector{T},
        compute_schedule,
        output_schedule,
        writer
    ) where T

Create a memory-efficient binned diagnostic using histogram binning.

# Arguments
- `variable`: The diagnostic variable to track distribution for
- `bin_edges`: Bin boundaries for the histogram (n+1 edges for n bins)
- `compute_schedule`: When to compute the diagnostic (typically every timestep)
- `output_schedule`: When to output the distribution statistics (e.g., daily, weekly, monthly)
- `writer`: Output writer for saving results

# Example: Monthly distribution of daily precipitation  
```julia
monthly_precip_dist = binned_diagnostic(
    precip_var,
    [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0],  # bins
    EveryDtSchedule(86400),                  # compute daily
    EveryDtSchedule(30*86400),               # output monthly
    NetCDFWriter(space, "output/")
)
```
"""
function binned_diagnostic(
    variable::DiagnosticVariable,
    bin_edges::Vector{T},
    compute_schedule,
    output_schedule,
    writer,
) where {T}

    # Validate inputs
    length(bin_edges) >= 2 || error("Need at least 2 bin edges")
    issorted(bin_edges) || error("Bin edges must be sorted")
    all(x -> x > 0, diff(bin_edges)) ||
        error("Bin edges must be strictly increasing")

    # Create binned-specific variable
    binned_var = DiagnosticVariable(
        compute = variable.compute,
        compute! = variable.compute!,
        short_name = "$(variable.short_name)_binned",
        long_name = "$(variable.long_name) binned",
        standard_name = variable.standard_name,
        units = variable.units,
        comments = "binned of $(variable.long_name) with $(length(bin_edges)-1) bins. $(variable.comments)",
    )

    # Create an instance of our callable reducer, which holds the initial state.
    initial_accumulator = BinnedAccumulator(bin_edges)
    reducer = BinnedReducer(initial_accumulator)

    # Create the scheduled diagnostic with binned reduction
    return ScheduledDiagnostic(
        variable = binned_var,
        compute_schedule_func = compute_schedule,
        output_schedule_func = output_schedule,
        reduction_time_func = reducer,
        pre_output_hook! = binned_pre_output_hook!,
        output_writer = ClimaDiagnostics.Writers.BinnedWriter(writer),
    )
end

"""
    binned_pre_output_hook!(accumulator_field, counter)

After writing output, this hook resets the accumulator for the next window.
It does so by creating new accumulators with zero counts, preserving bin edges.
"""
function binned_pre_output_hook!(accumulator_field, counter)
    @. accumulator_field =
        BinnedAccumulator(getproperty(accumulator_field, :bin_edges))
    return nothing
end
