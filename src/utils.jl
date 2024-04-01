function seconds_to_str_short(time::Real)
    name = ""
    days, rem_seconds = divrem(time, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    if round(seconds) == seconds
        seconds = convert(Int, seconds)
    end

    days > 0 && (name *= "$(days)d_")
    hours > 0 && (name *= "$(hours)h_")
    minutes > 0 && (name *= "$(minutes)m_")
    seconds > 0 && (name *= "$(seconds)s_")
    return rstrip(name, '_')
end

function seconds_to_str_long(time::Real)
    name = ""
    days, rem_seconds = divrem(time, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    days, hours, minutes = map(Int, (days, hours, minutes))

    if round(seconds) == seconds
        seconds = convert(Int, seconds)
    end

    days > 0 && (name *= "$(days) Day(s) ")
    hours > 0 && (name *= "$(hours) Hour(s) ")
    minutes > 0 && (name *= "$(minutes) Minute(s) ")
    seconds > 0 && (name *= "$(seconds) Second(s) ")

    return rstrip(name, ' ')
end
