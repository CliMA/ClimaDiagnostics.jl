module ScheduledDiagnostics

import ..AbstractWriter
import ..Callbacks: EveryStepSchedule
import ..DiagnosticVariables:
    DiagnosticVariable, descriptive_short_name, descriptive_long_name

"""
    ScheduledDiagnostic

Conceptually, a ScheduledDiagnostics is a DiagnosticVariable we want to compute in a given
simulation. For example, it could be the temperature averaged over a day. We can have
multiple ScheduledDiagnostics for the same DiagnosticVariable (e.g., daily and monthly
average temperatures).
"""
struct ScheduledDiagnostic{
    T1,
    T2,
    OW <: AbstractWriter,
    F1,
    PO,
    DV <: DiagnosticVariable,
}
    """The `DiagnosticVariable` that has to be computed and output"""
    variable::DV

    """A boolean function that determines when this diagnostic should be output. It has to
    take one argument, the integrator. Most typically, only `integrator.t` or
    `integrator.step` are used. Could be a Callback.AbstractSchedule."""
    output_schedule_func::T1

    """Struct that controls out to save the computed diagnostic variable to disk.
    `output_writer` has to implement a method `write_field!` that takes three arguments: the
    value that has to be output, the `ScheduledDiagnostic`, and the integrator. Internally,
    the integrator contains extra information (such as the current timestep). It is
    responsibility of the `output_writer` to properly use the provided information for
    meaningful output."""
    output_writer::OW

    """If not `nothing`, this `ScheduledDiagnostic` receives an area of scratch space `acc`
    where to accumulate partial results. Then, ar directed by the `compute_schedule_func`,
    `reduction_time_func` is computed between the previously stored value in `acc` and the
    new value. This implements a running reduction. For example, if `reduction_time_func =
    max`, the space `acc` will hold the running maxima of the diagnostic. To implement
    operations like the arithmetic average, the `reduction_time_func` has to be chosen as
    `sum`, and a `pre_output_hook!` that renormalizes `acc` by the number of samples has to
    be provided. For custom reductions, it is necessary to also specify the identity of
    operation by defining a new method to `identity_of_reduction`."""
    reduction_time_func::F1

    """A boolean function that determines when this diagnostic should be computed. It has to
    take one argument, the integrator. Most typically, only `integrator.t` or
    `integrator.step` are used. Could be a Callback.AbstractSchedule."""
    compute_schedule_func::T2

    # Design note: pre_output_hook!
    #
    # One of our key requirements is to be able to compute arithmetic averages.
    # Unfortunately, computing an arithmetic average requires keeping track of how many
    # elements we are summing up. pre_output_hook! was introduced so that we can think of an
    # average as a sum coupled with division, and perform the division (by the number of
    # elements) before output. pre_output_hook! could be used for other operations, but we
    # decided to keep it simple and target directly the most important use case for us.
    #
    # This choice restricts what reductions can be performed. For example, it is not
    # possible to have a geometric average. If more complex reduction are needed, this
    # mechanism has to be changed.

    """Function that has to be run before saving to disk for reductions (mostly used to
    implement averages). The function `pre_output_hook!` is called with two arguments: the value
    accumulated during the reduction, and the number of times the diagnostic was computed from
    the last time it was output. `pre_output_hook!` should mutate the accumulator in place. The
    return value of `pre_output_hook!` is discarded. An example of `pre_output_hook!` to compute
    the arithmetic average is `pre_output_hook!(acc, N) = @. acc = acc / N`."""
    pre_output_hook!::PO

    """Short name used to output this ScheduledDiagnostic for file names or datasets."""
    output_short_name::String

    """Descriptive name used to output this ScheduledDiagnostic in metadata or attributes."""
    output_long_name::String
end

"""
     ScheduledDiagnosticIterations(; variable::DiagnosticVariable,
                                     output_schedule_func,
                                     output_writer,
                                     reduction_time_func = nothing,
                                     reduction_space_func = nothing,
                                     compute_schedule_func = isa_reduction ? 1 : output_every,
                                     pre_output_hook! = nothing,
                                     output_short_name = descriptive_short_name(self),
                                     output_short_name = descriptive_long_name(self))


 A `DiagnosticVariable` that has to be computed and output during a simulation with a cadence
 defined by the number of iterations, with an optional reduction applied to it (e.g., compute
 the maximum temperature over the course of every 10 timesteps). This object is turned into
 two callbacks (one for computing and the other for output) and executed by the integrator.

 Keyword arguments
 =================

 - `variable`: The `DiagnosticVariable` that has to be computed and output.

 - `output_every`: Save the results to disk every `output_every` iterations. If `output_every`
                   is non-positive, only output at the first time step.

 - `output_writer`: Function that controls out to save the computed diagnostic variable to
                    disk. `output_writer` has to take three arguments: the value that has to
                    be output, the `ScheduledDiagnostic`, and the integrator. Internally, the
                    integrator contains extra information (such as the current timestep). It
                    is responsibility of the `output_writer` to properly use the provided
                    information for meaningful output.

 - `reduction_time_func`: If not `nothing`, this `ScheduledDiagnostic` receives an area of
                          scratch space `acc` where to accumulate partial results. Then, at
                          every `compute_every`, `reduction_time_func` is computed between
                          the previously stored value in `acc` and the new value. This
                          implements a running reduction. For example, if
                          `reduction_time_func = max`, the space `acc` will hold the running
                          maxima of the diagnostic. To implement operations like the
                          arithmetic average, the `reduction_time_func` has to be chosen as
                          `sum`, and a `pre_output_hook!` that renormalize `acc` by the
                          number of samples has to be provided. For custom reductions, it is
                          necessary to also specify the identity of operation by defining a
                          new method to `identity_of_reduction`.

 - `reduction_space_func`: NOT IMPLEMENTED YET

 - `compute_every`: Run the computations every `compute_every` iterations. This is not
                    particularly useful for point-wise diagnostics, where we enforce that
                    `compute_every` = `output_every`. For time reductions, `compute_every` is
                    set to 1 (compute at every timestep) by default. `compute_every` has to
                    evenly divide `output_every`.

 - `pre_output_hook!`: Function that has to be run before saving to disk for reductions
                       (mostly used to implement averages). The function `pre_output_hook!`
                       is called with two arguments: the value accumulated during the
                       reduction, and the number of times the diagnostic was computed from
                       the last time it was output. `pre_output_hook!` should mutate the
                       accumulator in place. The return value of `pre_output_hook!` is
                       discarded. An example of `pre_output_hook!` to compute the arithmetic
                       average is `pre_output_hook!(acc, N) = @. acc = acc / N`.

- `output_short_name`: A descriptive name for this particular diagnostic. If none is
                       provided, one will be generated mixing the short name of the
                       variable, the reduction, and the period of the reduction.
                       Normally, it has to be unique. In `ClimaAtmos`, we follow the CMIP
                       conventions for this.

- `output_long_name`: A descriptive name for this particular diagnostic. If none is
                      provided, one will be generated mixing the short name of the
                      variable, the reduction, and the period of the reduction.

 """
function ScheduledDiagnostic(;
    variable::DiagnosticVariable,
    output_writer,
    reduction_time_func = nothing,
    compute_schedule_func = EveryStepSchedule(),
    output_schedule_func = isnothing(reduction_time_func) ?
                           compute_schedule_func : EveryStepSchedule(),
    pre_output_hook! = (accum, count) -> nothing,
    output_short_name = descriptive_short_name(
        variable,
        output_schedule_func,
        reduction_time_func,
        pre_output_hook!,
    ),
    output_long_name = descriptive_long_name(
        variable,
        output_schedule_func,
        reduction_time_func,
        pre_output_hook!,
    ),
)
    # pre_output_hook! has to be a function, but it is much more intuitive to specify
    # `nothing` when we want nothing to happen. Here, we convert the nothing keyword
    # into a function that does nothing
    if isnothing(pre_output_hook!)
        pre_output_hook! = (accum, count) -> nothing
    end

    T = typeof(variable)
    T1 = typeof(output_schedule_func)
    T2 = typeof(compute_schedule_func)
    OW = typeof(output_writer)
    F1 = typeof(reduction_time_func)
    PO = typeof(pre_output_hook!)

    ScheduledDiagnostic{T1, T2, OW, F1, PO, T}(
        variable,
        output_schedule_func,
        output_writer,
        reduction_time_func,
        compute_schedule_func,
        pre_output_hook!,
        output_short_name,
        output_long_name,
    )
end

"""
    output_short_name(sd::ScheduledDiagnostic)

Return the short name to use for output of the `sd` `ScheduledDiagnostic`.
"""
function output_short_name(sd::ScheduledDiagnostic)
    return sd.output_short_name
end

# """
#     ScheduledDiagnosticIterations(; variable::DiagnosticVariable,
#                                     output_schedule_func,
#                                     output_writer,
#                                     reduction_time_func = nothing,
#                                     reduction_space_func = nothing,
#                                     compute_schedule_func = isa_reduction ? 1 : output_every,
#                                     pre_output_hook! = nothing,
#                                     output_short_name = descriptive_short_name(self),
#                                     output_short_name = descriptive_long_name(self))


# A `DiagnosticVariable` that has to be computed and output during a simulation with a cadence
# defined by the number of iterations, with an optional reduction applied to it (e.g., compute
# the maximum temperature over the course of every 10 timesteps). This object is turned into
# two callbacks (one for computing and the other for output) and executed by the integrator.

# Keyword arguments
# =================

# - `variable`: The diagnostic variable that has to be computed and output.

# - `output_every`: Save the results to disk every `output_every` iterations. If `output_every`
#                   is non-positive, only output at the first time step.

# - `output_writer`: Function that controls out to save the computed diagnostic variable to
#                    disk. `output_writer` has to take three arguments: the value that has to
#                    be output, the `ScheduledDiagnostic`, and the integrator. Internally, the
#                    integrator contains extra information (such as the current timestep). It
#                    is responsibility of the `output_writer` to properly use the provided
#                    information for meaningful output.

# - `reduction_time_func`: If not `nothing`, this `ScheduledDiagnostic` receives an area of
#                          scratch space `acc` where to accumulate partial results. Then, at
#                          every `compute_every`, `reduction_time_func` is computed between
#                          the previously stored value in `acc` and the new value. This
#                          implements a running reduction. For example, if
#                          `reduction_time_func = max`, the space `acc` will hold the running
#                          maxima of the diagnostic. To implement operations like the
#                          arithmetic average, the `reduction_time_func` has to be chosen as
#                          `sum`, and a `pre_output_hook!` that renormalize `acc` by the
#                          number of samples has to be provided. For custom reductions, it is
#                          necessary to also specify the identity of operation by defining a
#                          new method to `identity_of_reduction`.

# - `reduction_space_func`: NOT IMPLEMENTED YET

# - `compute_every`: Run the computations every `compute_every` iterations. This is not
#                    particularly useful for point-wise diagnostics, where we enforce that
#                    `compute_every` = `output_every`. For time reductions, `compute_every` is
#                    set to 1 (compute at every timestep) by default. `compute_every` has to
#                    evenly divide `output_every`.

# - `pre_output_hook!`: Function that has to be run before saving to disk for reductions
#                       (mostly used to implement averages). The function `pre_output_hook!`
#                       is called with two arguments: the value accumulated during the
#                       reduction, and the number of times the diagnostic was computed from
#                       the last time it was output. `pre_output_hook!` should mutate the
#                       accumulator in place. The return value of `pre_output_hook!` is
#                       discarded. An example of `pre_output_hook!` to compute the arithmetic
#                       average is `pre_output_hook!(acc, N) = @. acc = acc / N`.

#  `output_short_name`: A descriptive name for this particular diagnostic. If none is
#                       provided, one will be generated mixing the short name of the
#                       variable, the reduction, and the period of the reduction.
#                       Normally, it has to be unique. In `ClimaAtmos`, we follow the CMIP
#                       conventions for this.

#  `output_long_name`: A descriptive name for this particular diagnostic. If none is
#                      provided, one will be generated mixing the short name of the
#                      variable, the reduction, and the period of the reduction.

# """
# function ScheduledDiagnosticIterations(;
#     variable::DiagnosticVariable{T},
#     output_every,
#     output_writer,
#     reduction_time_func = nothing,
#     reduction_space_func = nothing,
#     compute_every = isnothing(reduction_time_func) ? output_every : 1,
#     pre_output_hook! = nothing,
#     output_short_name = descriptive_short_name(
#         variable,
#         output_every,
#         reduction_time_func,
#         pre_output_hook!;
#     ),
#     output_long_name = descriptive_long_name(
#         variable,
#         output_every,
#         reduction_time_func,
#         pre_output_hook!;
#     ),
# ) where {T}

#     # We provide an inner constructor to enforce some constraints

#     (output_every <= 0 || output_every % compute_every == 0) || error(
#         "output_every ($output_every) should be multiple of compute_every ($compute_every) for diagnostic $(output_short_name)",
#     )

#     isa_reduction = !isnothing(reduction_time_func)

#     # If it is not a reduction, we compute only when we output
#     if !isa_reduction && compute_every != output_every
#         @warn "output_every ($output_every) != compute_every ($compute_every) for $(output_short_name), changing compute_every to match"
#         compute_every = output_every
#     end

#     # pre_output_hook! has to be a function, but it is much more intuitive to specify
#     # `nothing` when we want nothing to happen. Here, we convert the nothing keyword
#     # into a function that does nothing
#     if isnothing(pre_output_hook!)
#         pre_output_hook! = (accum, count) -> nothing
#     end

#     T1 = typeof(output_every)
#     T2 = typeof(compute_every)
#     OW = typeof(output_writer)
#     F1 = typeof(reduction_time_func)
#     F2 = typeof(reduction_space_func)
#     PO = typeof(pre_output_hook!)

#     new{T1, T2, OW, F1, F2, PO, T}(
#         variable,
#         output_every,
#         output_writer,
#         reduction_time_func,
#         reduction_space_func,
#         compute_every,
#         pre_output_hook!,
#         output_short_name,
#         output_long_name,
#     )
# end



# struct ScheduledDiagnosticTime{T1, T2, OW, F1, F2, PO}
#     variable::DiagnosticVariable
#     output_every::T1
#     output_writer::OW
#     reduction_time_func::F1
#     reduction_space_func::F2
#     compute_every::T2
#     pre_output_hook!::PO
#     output_short_name::String
#     output_long_name::String

#     """
#         ScheduledDiagnosticTime(; variable::DiagnosticVariable,
#                                   output_every,
#                                   output_writer,
#                                   reduction_time_func = nothing,
#                                   reduction_space_func = nothing,
#                                   compute_every = isa_reduction ? :timestep : output_every,
#                                   pre_output_hook! = nothing,
#                                   output_short_name = descriptive_short_name(self),
#                                   output_long_name = descriptive_long_name(self),
#                                   )


#     A `DiagnosticVariable` that has to be computed and output during a simulation with a
#     cadence defined by how many seconds in simulation time, with an optional reduction
#     applied to it (e.g., compute the maximum temperature over the course of every day). This
#     object is turned into a `ScheduledDiagnosticIterations`, which is turned into two
#     callbacks (one for computing and the other for output) and executed by the integrator.

#     Keyword arguments
#     =================

#     - `variable`: The diagnostic variable that has to be computed and output.

#     - `output_every`: Save the results to disk every `output_every` seconds. If `output_every`
#                       is non-positive, only output at the first time step.

#     - `output_writer`: Function that controls out to save the computed diagnostic variable to
#                        disk. `output_writer` has to take three arguments: the value that has to
#                        be output, the `ScheduledDiagnostic`, and the integrator. Internally, the
#                        integrator contains extra information (such as the current timestep). It
#                        is responsibility of the `output_writer` to properly use the provided
#                        information for meaningful output.

#     - `reduction_time_func`: If not `nothing`, this `ScheduledDiagnostic` receives an area of
#                              scratch space `acc` where to accumulate partial results. Then, at
#                              every `compute_every`, `reduction_time_func` is computed between
#                              the previously stored value in `acc` and the new value. This
#                              implements a running reduction. For example, if
#                              `reduction_time_func = max`, the space `acc` will hold the running
#                              maxima of the diagnostic. To implement operations like the
#                              arithmetic average, the `reduction_time_func` has to be chosen as
#                              `sum`, and a `pre_output_hook!` that renormalize `acc` by the
#                              number of samples has to be provided. For custom reductions, it is
#                              necessary to also specify the identity of operation by defining a
#                              new method to `identity_of_reduction`.

#     - `reduction_space_func`: NOT IMPLEMENTED YET

#     - `compute_every`: Run the computations every `compute_every` seconds. This is not
#                        particularly useful for point-wise diagnostics, where we enforce that
#                        `compute_every` = `output_every`. For time reductions,
#                        `compute_every` is set to `:timestep` (compute at every timestep) by
#                        default. `compute_every` has to evenly divide `output_every`.
#                        `compute_every` can take the special symbol `:timestep` which is a
#                        placeholder for the timestep of the simulation to which this
#                        `ScheduledDiagnostic` is attached.

#     - `pre_output_hook!`: Function that has to be run before saving to disk for reductions
#                           (mostly used to implement averages). The function `pre_output_hook!`
#                           is called with two arguments: the value accumulated during the
#                           reduction, and the number of times the diagnostic was computed from
#                           the last time it was output. `pre_output_hook!` should mutate the
#                           accumulator in place. The return value of `pre_output_hook!` is
#                           discarded. An example of `pre_output_hook!` to compute the arithmetic
#                           average is `pre_output_hook!(acc, N) = @. acc = acc / N`.

#    - `output_short_name`: A descriptive name for this particular diagnostic. If none is
#                           provided, one will be generated mixing the short name of the
#                           variable, the reduction, and the period of the reduction.
#                           Normally, it has to be unique. In `ClimaAtmos`, we follow the CMIP
#                           conventions for this.

#    - `output_long_name`: A descriptive name for this particular diagnostic. If none is
#                          provided, one will be generated mixing the short name of the
#                          variable, the reduction, and the period of the reduction.
#     """
#     function ScheduledDiagnosticTime(;
#         variable::DiagnosticVariable,
#         output_every,
#         output_writer,
#         reduction_time_func = nothing,
#         reduction_space_func = nothing,
#         compute_every = isnothing(reduction_time_func) ? output_every :
#                         :timestep,
#         pre_output_hook! = nothing,
#         output_short_name = descriptive_short_name(
#             variable,
#             output_every,
#             reduction_time_func,
#             pre_output_hook!;
#         ),
#         output_long_name = descriptive_long_name(
#             variable,
#             output_every,
#             reduction_time_func,
#             pre_output_hook!;
#         ),
#     )

#         # We provide an inner constructor to enforce some constraints

#         # compute_every could be a Symbol (:timestep). We process this that when we process
#         # the list of diagnostics
#         if !isa(compute_every, Symbol)
#             (output_every <= 0 || output_every % compute_every == 0) || error(
#                 "output_every ($output_every) should be multiple of compute_every ($compute_every) for diagnostic $(output_short_name)",
#             )
#         end

#         isa_reduction = !isnothing(reduction_time_func)

#         # If it is not a reduction, we compute only when we output
#         if !isa_reduction && compute_every != output_every
#             @warn "output_every ($output_every) != compute_every ($compute_every) for $(output_short_name), changing compute_every to match"
#             compute_every = output_every
#         end

#         # pre_output_hook! has to be a function, but it is much more intuitive to specify
#         # `nothing` when we want nothing to happen. Here, we convert the nothing keyword
#         # into a function that does nothing
#         if isnothing(pre_output_hook!)
#             pre_output_hook! = (accum, count) -> nothing
#         end

#         T1 = typeof(output_every)
#         T2 = typeof(compute_every)
#         OW = typeof(output_writer)
#         F1 = typeof(reduction_time_func)
#         F2 = typeof(reduction_space_func)
#         PO = typeof(pre_output_hook!)

#         new{T1, T2, OW, F1, F2, PO}(
#             variable,
#             output_every,
#             output_writer,
#             reduction_time_func,
#             reduction_space_func,
#             compute_every,
#             pre_output_hook!,
#             output_short_name,
#             output_long_name,
#         )
#     end
# end

# """
#     ScheduledDiagnosticIterations(sd_time::ScheduledDiagnosticTime, Δt)


# Create a `ScheduledDiagnosticIterations` given a `ScheduledDiagnosticTime` and a timestep
# `Δt`. In this, ensure that `compute_every` and `output_every` are meaningful for the given
# timestep.

# """
# function ScheduledDiagnosticIterations(
#     sd_time::ScheduledDiagnosticTime,
#     Δt::T,
# ) where {T}

#     # If we have the timestep, we can convert time in seconds into iterations

#     # if compute_every is :timestep, then we want to compute after every iterations
#     compute_every =
#         sd_time.compute_every == :timestep ? 1 : sd_time.compute_every / Δt
#     output_every = sd_time.output_every / Δt

#     # When Δt is a Float32, loss of precision might lead to spurious results (e.g., 1. /
#     # 0.1f0 = 9.99999985098839). So, we round to the number of significant digits that we
#     # expect from the float type.
#     #
#     # FIXME: eps(typeof(Δt)) is not the best value to pick the number of significant digits
#     # because it makes sense only for values of order unity.
#     sigdigits = eps(typeof(Δt)) |> log10 |> abs |> round |> Int

#     output_every = round(output_every; sigdigits)
#     compute_every = round(compute_every; sigdigits)

#     isinteger(output_every) || error(
#         "output_every ($(sd_time.output_every)) should be multiple of the timestep ($Δt) for diagnostic $(sd_time.output_short_name)",
#     )
#     isinteger(compute_every) || error(
#         "compute_every ($(sd_time.compute_every)) should be multiple of the timestep ($Δt) for diagnostic $(sd_time.output_short_name)",
#     )

#     ScheduledDiagnosticIterations(;
#         sd_time.variable,
#         output_every = convert(Int, output_every),
#         sd_time.output_writer,
#         sd_time.reduction_time_func,
#         sd_time.reduction_space_func,
#         compute_every = convert(Int, compute_every),
#         sd_time.pre_output_hook!,
#         sd_time.output_short_name,
#         sd_time.output_long_name,
#     )
# end

# """
#     ScheduledDiagnosticTime(sd_time::ScheduledDiagnosticIterations, Δt)


# Create a `ScheduledDiagnosticTime` given a `ScheduledDiagnosticIterations` and a timestep
# `Δt`.

# """
# function ScheduledDiagnosticTime(
#     sd_time::ScheduledDiagnosticIterations,
#     Δt::T,
# ) where {T}

#     # If we have the timestep, we can convert time in iterations to seconds

#     # if compute_every is :timestep, then we want to compute after every iterations
#     compute_every =
#         sd_time.compute_every == 1 ? :timestep : sd_time.compute_every * Δt
#     output_every = sd_time.output_every * Δt

#     ScheduledDiagnosticTime(;
#         sd_time.variable,
#         output_every,
#         sd_time.output_writer,
#         sd_time.reduction_time_func,
#         sd_time.reduction_space_func,
#         compute_every,
#         sd_time.pre_output_hook!,
#         sd_time.output_short_name,
#         sd_time.output_long_name,
#     )
# end

# # We provide also a companion constructor for ScheduledDiagnosticIterations which returns
# # itself (without copy) when called with a timestep.
# #
# # This is so that we can assume that
# # ScheduledDiagnosticIterations(ScheduledDiagnostic{Time, Iterations}, Δt)
# # always returns a valid ScheduledDiagnosticIterations
# ScheduledDiagnosticIterations(
#     sd::ScheduledDiagnosticIterations,
#     _Δt::T,
# ) where {T} = sd
end
