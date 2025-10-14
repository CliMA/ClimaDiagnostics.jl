"""
    NetCDFBuffer{
        NC <: NCDatasets.NCDataset,
        DATE <: Union{Nothing, Dates.AbstractDateTime},
        TT,
        ARRAY <: Array,
    }

Store times and data that will eventually be written to the corresponding NetCDF file.

This struct is not responsible for populating the NetCDF file with the appropriate
dimensions and variable. This struct is responsible for caching and writing values of the
time, date, time bounds, and date bounds dimension, and the data of the variable.
"""
struct NetCDFBuffer{
    NC <: NCDatasets.NCDataset,
    DATE <: Union{Nothing, Dates.AbstractDateTime},
    TT,
    ARRAY <: Array,
}
    """NetCDF file to write to when flushing data in buffer"""
    nc::NC

    """Short name of the variable"""
    short_name::String

    """Start date of the simulation"""
    start_date::DATE

    """A reference to the index along the time dimension that the data must be written to in
    the NetCDF file"""
    t_start_idx::Base.RefValue{Int64}

    """Store all the times, either floats or `ITime`s, before writing to NetCDF file."""
    times::Vector{TT}

    """A buffer for storing data."""
    storage::ARRAY

    """The max capacity of the buffer."""
    max_capacity::Int64
end

"""
    NetCDFBuffer(
        t,
        data;
        max_capacity = 1,
    )

Construct a `NetCDFBuffer` of size `max_capacity`.

The buffer is initialized with time `t` in the times and `data` in the storage.
"""
function NetCDFBuffer(
    nc::NCDatasets.NCDataset,
    short_name::String,
    start_date::Union{Nothing, Dates.AbstractDateTime},
    t,
    data;
    max_capacity::Integer = 1,
)
    max_capacity < 1 && error(
        "Cannot construct NetCDFBuffer with a max_capacity of $max_capacity. Provide a positive number for max_capacity.",
    )
    # When flushing the data, the time index to write the data to begin at 1
    t_start_idx = Base.RefValue(1)
    times = [t]
    sizehint!(times, max_capacity)

    storage_size = (max_capacity, size(data)...)
    storage = similar(data, storage_size)
    selectdim(storage, 1, 1) .= data
    return NetCDFBuffer{
        typeof(nc),
        typeof(start_date),
        eltype(times),
        typeof(storage),
    }(
        nc,
        short_name,
        start_date,
        t_start_idx,
        times,
        storage,
        max_capacity,
    )
end

"""
    empty!(buffer::NetCDFBuffer)

Empty the `buffer` and return nothing.

This does not write to the NetCDF file.
"""
function Base.empty!(buffer::NetCDFBuffer)
    buffer.t_start_idx[] += length(buffer.times)
    empty!(buffer.times)
    return nothing
end

"""
    isempty(buffer::NetCDFBuffer)

Return whether the `buffer` is empty.
"""
function Base.isempty(buffer::NetCDFBuffer)
    return isempty(buffer.times)
end

"""
    length(buffer::NetCDFBuffer)

Return the number of `times` in the buffer.
"""
function Base.length(buffer::NetCDFBuffer)
    return length(buffer.times)
end

"""
    isfull(buffer::NetCDFBuffer)

Return whether the `buffer` is full.
"""
function isfull(buffer::NetCDFBuffer)
    return length(buffer) == buffer.max_capacity
end

"""
    push!(buffer::NetCDFBuffer, time, data)

Push `time` and `data` to the buffer.

If the `buffer` is full, then `flush!` is called before writing to the `buffer`.

As a result of flushing when the buffer is full, it is not recommended to use a buffer of
size 1.
"""
function Base.push!(buffer::NetCDFBuffer, time, data)
    isfull(buffer) && flush!(buffer)
    (; times, storage) = buffer
    push!(times, time)
    selectdim(storage, 1, length(times)) .= data
    return nothing
end

"""
    flush!(buffer::NetCDFBuffer)

Flush the buffer by writing the times, dates, and data to the NetCDF file and emptying the
`buffer`.
"""
function flush!(buffer::NetCDFBuffer)
    (; nc, short_name, start_date, t_start_idx, times, storage) = buffer
    length(buffer) >= 1 || return nothing
    seconds = float.(times)

    # Compute time indices from t_start_idx
    t_start_idx = t_start_idx[]
    time_range = t_start_idx:(t_start_idx + length(seconds) - 1)

    _write_temporal_dims!(nc, "time", "time_bnds", seconds, time_range)

    if !isnothing(start_date)
        dates = start_date .+ Dates.Millisecond.(round.(1000 * seconds))
        _write_temporal_dims!(
            nc,
            "date",
            "date_bnds",
            dates,
            time_range,
            start_date = start_date,
        )
    end

    colons = ntuple(_ -> Colon(), ndims(storage) - 1)
    interpolate_field = selectdim(storage, 1, eachindex(seconds))
    nc[short_name][time_range, colons...] = interpolate_field
    empty!(buffer)
    return nothing
end

"""
    _write_temporal_dims!(
        nc,
        time_name,
        time_bnds_name,
        times,
        time_range;
        start_date = nothing,
    )

Write the temporal dimensions for either the `time` and `time_bnds` dimensions or `date` and
`date_bnds` dimensions.
"""
function _write_temporal_dims!(
    nc,
    time_name,
    time_bnds_name,
    times,
    time_range;
    start_date = nothing,
)
    first_t_idx = first(time_range)
    first_time = first(times)
    nc[time_name][time_range] = times
    if first_t_idx == 1
        first_lower_time_bnd = ifelse(
            first(times) isa Dates.Dates.AbstractDateTime,
            start_date,
            zero(first_time),
        )
        nc[time_bnds_name][:, 1] = [first_lower_time_bnd; first_time]
    else
        nc[time_bnds_name][:, first_t_idx] =
            [nc[time_name][first_t_idx - 1]; first_time]
    end
    if length(times) > 1
        second_t_idx = first_t_idx + 1
        last_t_idx = last(time_range)
        nc[time_bnds_name][1, second_t_idx:last_t_idx] = times[1:(end - 1)]
        nc[time_bnds_name][2, second_t_idx:last_t_idx] = times[2:end]
    end
    return nothing
end
