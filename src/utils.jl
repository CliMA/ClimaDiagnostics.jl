"""
    seconds_to_str_short(time_seconds::Real)

Convert a time in seconds to a string representing the time in "short" units.

Examples:
========

```julia
julia> seconds_to_str_short(0)
"0s"

julia> seconds_to_str_short(0.1)
"0.1s"

julia> seconds_to_str_short(60)
"1m"

julia> seconds_to_str_short(3600)
"1h"

julia> seconds_to_str_short(3601)
"1h_1s"
```
"""
function seconds_to_str_short(time_seconds::Real)
    name = ""
    days, rem_seconds = divrem(time_seconds, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    round(seconds) ≈ seconds && (seconds = convert(Int, seconds))

    days > 0 && (name *= "$(days)d_")
    hours > 0 && (name *= "$(hours)h_")
    minutes > 0 && (name *= "$(minutes)m_")
    seconds > 0 && (name *= "$(seconds)s_")
    return rstrip(name, '_')
end

"""
    seconds_to_str_long(time_seconds::Real)

Convert a time in seconds to a string representing the time in "longer" units.

Examples:
========

```julia
julia> seconds_to_str_long(0)
"0 Seconds"

julia> seconds_to_str_long(60)
"1 Minutes"

julia> seconds_to_str_long(3600)
"1 Hours"

julia> seconds_to_str_long(3601)
"1 Hours 1 Seconds"
```
"""
function seconds_to_str_long(time_seconds::Real)
    name = ""
    days, rem_seconds = divrem(time_seconds, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    round(seconds) ≈ seconds && (seconds = convert(Int, seconds))

    days > 0 && (name *= "$(days) Days ")
    hours > 0 && (name *= "$(hours) Hours ")
    minutes > 0 && (name *= "$(minutes) Minutes ")
    seconds > 0 && (name *= "$(seconds) Seconds ")

    return rstrip(name, ' ')
end
