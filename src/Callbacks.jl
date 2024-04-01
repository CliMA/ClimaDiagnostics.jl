# Callback.jl
module Callbacks

"""
    Callback

A lightweight struct that contains two functions, `callback_func`, the function that has to
be called, and `schedule_func`, a boolean function that determines if it is time to call
`func` or not.
"""
struct Callback{FUNC <: Function, SCHEDULE}
    """Function to be called. It has to take one argument, the integrator."""
    callback_func::FUNC

    """Boolean function (or, more often, a callable struct) that determines whether
    `callback_func` should be called or not. It has to take one argument, the integrator.
    Most typically, only `integrator.t` or `integrator.step` are used."""
    schedule_func::SCHEDULE
end

#############
# Schedules #
#############

"""
    DivisorSchedule

True when the iteration number is evenly divisible by a given number (this is roughly
equivalent to: "run this call back every N steps").
"""
struct DivisorSchedule
    """Return true when the step number is divided evenly by this number (ie, step % divisor
        == 0) """
    divisor::Int

    """String that can be used to identify this schedule when saving files/datasets."""
    filename_str::String

    function DivisorSchedule(divisor::Int)
        filename_str = string(divisor)
        new(divisor, filename_str)
    end
end

function (schedule::DivisorSchedule)(integrator)
    return rem(integrator.step, schedule.divisor) == 0
end

end
