module DiagnosticVariables

import ..Schedules: AbstractSchedule, long_name

"""
    DiagnosticVariable(;
        compute!,
        short_name = "",
        long_name = "",
        standard_name = "",
        units = "",
        comments = ""
    )

A recipe to compute a diagnostic variable from the state, along with some useful metadata.

The primary use for `DiagnosticVariable`s is to be embedded in a `ScheduledDiagnostic` to
compute diagnostics while the simulation is running.

The metadata is used exclusively by the `output_writer` in the `ScheduledDiagnostic`. It is
responsibility of the `output_writer` to follow the conventions about the meaning of the
metadata and their use.

In `ClimaAtmos`, we roughly follow the naming conventions listed in this file:
https://airtable.com/appYNLuWqAgzLbhSq/shrKcLEdssxb8Yvcp/tblL7dJkC3vl5zQLb

!!! compat "ClimaDiagnostics 0.2.13"

    Support for `compute` was introduced in version `0.2.13`. Prior to this version, the
    in-place `compute!` had to be provided.

Keyword arguments
=================

- `compute!`: Function that computes the diagnostic variable from the `state`, `cache`, and
              `time`. In addition to these three arguments, `compute!` has to take four
              arguments, the destination where to write the result of the computation.

- `compute`: Function that computes the diagnostic variable from the `state`, `cache`, and
             `time` and returns either a `Field` or an un-evaluated expression (e.g., with
             `LazyBroadcast.lazy`).

- `short_name`: Name used to identify the variable in the output files and in the file
                names. Short but descriptive. `ClimaAtmos` follows the CMIP conventions and
                the diagnostics are identified by the short name.

- `long_name`: Name used to describe the variable in the output files.

- `standard_name`: Standard name, as in
                   http://cfconventions.org/Data/cf-standard-names/71/build/cf-standard-name-table.html

- `units`: Physical units of the variable. For a diagnostic whose `compute`/`compute!`
           returns a `Field` of `NamedTuple`s (one scalar diagnostic per entry), `units`
           may also be a `NamedTuple` mapping each component to its own units string, or a
           function `property_chain -> units` (where `property_chain` is a tuple of
           `Symbol`s identifying the component, e.g. `(:warm, :autoconversion)`). A plain
           string is applied to every component.

- `comments`: More verbose explanation of what the variable is, or comments related to how
              it is defined or computed.
"""
struct DiagnosticVariable{
    F1 <: Union{Function, Nothing},
    F2 <: Union{Function, Nothing},
}
    compute!::F1
    compute::F2
    short_name::String
    long_name::String
    standard_name::String
    units::Union{AbstractString, NamedTuple, Function}
    comments::String
end

function DiagnosticVariable(;
    compute = nothing,
    compute! = nothing,
    short_name = "",
    long_name = "",
    standard_name = "",
    units = "",
    comments = "",
)
    only_one_compute_provided = xor(isnothing(compute), isnothing(compute!))

    if !only_one_compute_provided
        error(
            "One and only one between `compute` and `compute!` has to be provided",
        )
    end

    return DiagnosticVariable(
        compute!,
        compute,
        short_name,
        long_name,
        standard_name,
        units,
        comments,
    )
end

"""
    short_name(dv::DiagnosticVariable)

Return the short name associated to the given `DiagnosticVariable`.
"""
function short_name(dv::DiagnosticVariable)
    return dv.short_name
end

"""
    long_name(dv::DiagnosticVariable)

Return the long name associated to the given `DiagnosticVariable`.
"""
function long_name(dv::DiagnosticVariable)
    return dv.long_name
end

"""
    average_pre_output_hook!

Function to use as `pre_output_hook!` for a `ScheduledDiagnostic` to compute an arithmetic average.
"""
function average_pre_output_hook!(accum, counter)
    @. accum = accum / counter
    return nothing
end

"""
    descriptive_short_name(variable::DiagnosticVariable,
                           output_schedule_func,
                           reduction_time_func,
                           pre_output_hook!)

Return a compact, unique-ish, identifier generated from the given information. This function
is useful for filenames and error messages.
"""
function descriptive_short_name(
    variable::DiagnosticVariable,
    output_schedule_func,
    reduction_time_func,
    pre_output_hook!;
)
    var = "$(variable.short_name)"
    isa_reduction = !isnothing(reduction_time_func)

    if isa_reduction
        red = "$(reduction_time_func)"

        # Let's check if we are computing the average. Note that this might slip under the
        # radar if the user passes their own pre_output_hook!.
        if reduction_time_func == (+) &&
           nameof(pre_output_hook!) == :average_pre_output_hook!
            red = "average"
        end
        suffix = "$red"
    else
        suffix = "inst"
    end
    return "$(var)_$(output_schedule_func)_$(suffix)"
end

"""
    descriptive_long_name(variable::DiagnosticVariable,
                          output_every,
                          reduction_time_func,
                          pre_output_hook!)

Return a verbose description of the given output variable.

This function is useful for attributes in output files.
"""
function descriptive_long_name(
    variable::DiagnosticVariable,
    output_schedule_func,
    reduction_time_func,
    pre_output_hook!;
)
    var = "$(variable.long_name)"
    isa_reduction = !isnothing(reduction_time_func)

    if isa_reduction
        red = "$(reduction_time_func)"

        # Let's check if we are computing the average. Note that this might slip under the
        # radar if the user passes their own pre_output_hook!.
        if reduction_time_func == (+) &&
           pre_output_hook! == average_pre_output_hook!
            red = "average"
        end

        if output_schedule_func isa AbstractSchedule
            suffix = "$(red) within $(long_name(output_schedule_func))"
        else
            suffix = red
        end
    else
        suffix = "Instantaneous"
    end
    return "$(var), $(suffix)"
end

"""
    component_units(units, property_chain)

Resolve the units for one scalar component of a (possibly composite,
i.e. `NamedTuple`-valued) diagnostic field.

`property_chain` is the `Tuple` of accessors identifying the component, as
produced by `ClimaCore.Fields.property_chains`: `Symbol`s for `NamedTuple`
field names and `Int`s for `Tuple`/`NTuple` element indices. Its length is
the nesting depth, and it is the empty tuple `()` for a plain scalar field.

- a `String` (or any `AbstractString`) is used for every component;
- a `Function` is called as `units(property_chain)`; it should return the
  units string, or `nothing` if the component has no units;
- a `NamedTuple` is descended following `property_chain` (supporting nested
  `NamedTuple`s); as a convenience, a flat `NamedTuple` keyed by the leaf
  name is also accepted.

Returns `nothing` when there is no entry for the component (callers turn this
into an empty units field; [`units_warnings`](@ref) reports it). An entry
explicitly set to `""` is returned as `""`: a deliberate empty value, not a
mismatch.
"""
component_units(units::AbstractString, _) = units
component_units(units, property_chain) = units(property_chain)
function component_units(units::NamedTuple, property_chain)
    isempty(property_chain) && return units
    key = first(property_chain)
    rest = Base.tail(property_chain)
    # `haskey(::NamedTuple, ::Int)` is positional indexing — an `Int` chain
    # element (a `Tuple`/`NTuple` index) must NOT match a `NamedTuple` key.
    # Guard both branches on `Symbol`; a non-`Symbol` resolves to `nothing`.
    if key isa Symbol && haskey(units, key)
        sub = getfield(units, key)
        return isempty(rest) ? sub : component_units(sub, rest)
    end
    leaf = last(property_chain)
    if leaf isa Symbol && haskey(units, leaf)
        return getfield(units, leaf)
    end
    return nothing
end

"""
    units_warnings(units, property_chains; name = "")

Emit a one-time warning listing the components in `property_chains` whose
`units` resolve to `nothing` (no entry) - their units field is left empty.

A `String` never resolves to `nothing` (it is shared by every component, and
the default `units` is `""`), so it is implicitly exempt. For a `Function`,
returning `nothing` signals "no units" (warned) while returning `""` is a
deliberate empty value (not warned); likewise a `NamedTuple` entry explicitly
set to `""` is deliberate and not warned about.

Intended to be called once per diagnostic (at `DiagnosticsHandler`
construction).
"""
function units_warnings(units, property_chains; name = "")
    no_units =
        filter(c -> isnothing(component_units(units, c)), property_chains)
    isempty(no_units) && return nothing
    label = isempty(name) ? "" : " for \"$name\""
    components = join((join(string.(c), '.') for c in no_units), ", ")
    @warn "Missing `units`$label for component(s): $components. \
           Their units attribute will be left empty."
    return nothing
end

"""
    units_attribute(units, property_chains)

Render `units` as a single string, for a writer that stores a (possibly
`NamedTuple`-valued) field as one dataset and therefore has a single `units`
attribute for the whole field (e.g. the `HDF5Writer`).

A plain `String` is returned unchanged. Otherwise the units of every component
in `property_chains` are resolved with [`component_units`](@ref) and joined,
e.g. `"warm.au: 1/s, warm.ac: 1/s, tot: 1/s"`. This keeps the attribute
consistent with the per-variable `units` the `NetCDFWriter` writes.
"""
units_attribute(units::AbstractString, _) = units
function units_attribute(units, property_chains)
    labeled = map(property_chains) do chain
        u = something(component_units(units, chain), "")
        "$(join(string.(chain), '.')): $u"
    end
    return join(labeled, ", ")
end
end
