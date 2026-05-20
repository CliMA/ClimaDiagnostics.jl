# NamedTuple diagnostics (one variable per component)

Sometimes a single computation naturally produces *several* related scalar
fields at once. A typical example is the individual process-level tendencies of a
parameterization as a `NamedTuple` at every grid point (e.g. condensation,
evaporation, accretion, ...).

`ClimaDiagnostics` lets you expose such a computation as a single
`DiagnosticVariable` whose `compute`/`compute!` returns a `Field` whose
*pointwise element type is a `NamedTuple`*. Each entry of the `NamedTuple` is
then written automatically as its own scalar diagnostic — you do not need to
write one `compute` function per component.

```julia
using ClimaDiagnostics

# `cache.microphysics_tendencies` is a Field of NamedTuples, e.g. with element
# type @NamedTuple{cond::FT, evap::FT, accr::FT}
compute_micro(state, cache, time) = cache.microphysics_tendencies

var = DiagnosticVariable(;
    # `short_name`/`long_name` are a *stem*, not the name of a single
    # variable: one diagnostic is produced per `NamedTuple` entry.
    short_name = "clw_tend",
    long_name = "cloud liquid water tendency",
    units = "kg kg⁻¹ s⁻¹",
    compute = compute_micro,
)
```

This single `var` does **not** describe one variable called `clw_tend`. It is
a recipe that fans out into one diagnostic per `NamedTuple` entry. With a
`NetCDFWriter` you get three variables,

| variable        | `long_name`                          | `units`       |
|-----------------|--------------------------------------|---------------|
| `clw_tend_cond` | `cloud liquid water tendency (cond)` | `kg kg⁻¹ s⁻¹` |
| `clw_tend_evap` | `cloud liquid water tendency (evap)` | `kg kg⁻¹ s⁻¹` |
| `clw_tend_accr` | `cloud liquid water tendency (accr)` | `kg kg⁻¹ s⁻¹` |

The per-component part of each name (`cond`, `evap`, `accr`) comes from the
keys of the `NamedTuple` your `compute` returns — you control it there, not
through `short_name`. `short_name`/`long_name` stay single strings: the stem
groups the components under a common prefix and (for the `NetCDFWriter`) also
names the single output file.

Wrapping `var` in a `ScheduledDiagnostic` and attaching a `NetCDFWriter` is
exactly the same as for an ordinary scalar diagnostic — see the [user
guide](@ref user_guide_header). Nothing else in the pipeline changes:
reductions, time averaging, schedules, and `ScheduledDiagnostic`s all work
unchanged (the accumulation operates element-wise on the `NamedTuple`).

## What gets written

The behavior depends on the writer, because only the `NetCDFWriter` needs to
flatten the field (see [Why splitting happens at the writer](@ref)).

### `NetCDFWriter`

The diagnostic is split into **one scalar NetCDF variable per component**, all
written to the same file. Each variable is named
`<short_name>_<component_path>`. With the `var` above you get the variables

```
clw_tend_cond
clw_tend_evap
clw_tend_accr
```

each interpolated and written exactly like an ordinary scalar diagnostic
(same dimensions, `time`/`date` bounds, reductions, ...). A plain scalar
diagnostic is unaffected: it keeps its bare `short_name` with no suffix.

### `HDF5Writer` and `DictWriter`

These writers do not interpolate, so they store the `NamedTuple` field
**whole and unsplit** — `DictWriter` keeps the `Field` as-is, and `HDF5Writer`
writes it as a single dataset that can be read back with `ClimaCore`'s
`InputOutput` module (the same mechanism used for the model state). No
per-component splitting happens.

## Per-component metadata

Only `units` can be specified per component. It accepts three forms:

- a `String` is applied to **every** component (the default, unchanged behavior);
- a `NamedTuple` maps each component to its own units string (it may mirror
  the nesting; a flat `NamedTuple` keyed by the leaf name also works);
- a function `property_chain -> units` is called for each component.
  `property_chain` is the `Tuple` identifying the component (`Symbol`s for
  `NamedTuple` field names, `Int`s for tuple-element indices, e.g.
  `(:warm, :autoconversion)`); the function should return the units string,
  or `nothing` if that component has no units. This is the most general form
  and is handy for dynamically generated or deeply nested results.

```julia
# All three forms are valid:
units = "kg kg⁻¹ s⁻¹"
units = (; cond = "kg kg⁻¹ s⁻¹", evap = "kg kg⁻¹ s⁻¹")
units = chain -> "kg kg⁻¹ s⁻¹"
```

`short_name` and `long_name` are deliberately single strings, **not**
per-component: `short_name` is a stem that groups the components and names the
output file, and the per-component identity is already carried by the
`NamedTuple` keys (the property chain). The `long_name` of each NetCDF
variable is the diagnostic's `long_name` with the component path appended,
e.g. `"cloud liquid water tendency (cond)"`.

For `HDF5Writer` (which keeps the field whole) the per-component units are
flattened into a single descriptive `variable_units` attribute, e.g.
`"cond: kg kg⁻¹ s⁻¹, evap: kg kg⁻¹ s⁻¹, accr: kg kg⁻¹ s⁻¹"`.

### Misspecified units

If a component has no units, its units field is left empty (an empty NetCDF
`units` attribute / an empty value in the `HDF5Writer`'s flattened string),
and a one-time warning listing those components is emitted when the
`DiagnosticsHandler` is constructed.

A component "has no units" when it resolves to `nothing`: a `NamedTuple` with
no entry for it, or a `Function` that returned `nothing` for it. An
explicitly empty value is treated as deliberate and is **not** warned about:
a `NamedTuple` entry set to `""`, or a `Function` returning `""`. A `String`
is never warned (it is shared by every component, and the default `units` is
`""`).

## Nested `NamedTuple`s

Nesting is supported and flattened recursively. A field with element type

```julia
@NamedTuple{warm::@NamedTuple{au::FT, ac::FT}, tot::FT}
```

and `short_name = "proc"` produces the NetCDF variables

```
proc_warm_au
proc_warm_ac
proc_tot
```

The corresponding property chains are `(:warm, :au)`, `(:warm, :ac)`, and
`(:tot,)`; these are exactly the chains passed to a `units` function and used
to build the flat-`NamedTuple` lookup.

## Point and column spaces

Single-column and `PointSpace` diagnostics are supported too: the
`NetCDFWriter` splits the `NamedTuple` into one scalar variable per component
on these spaces as well (no interpolation is performed).

## Reductions and averaging

Temporal reductions work without any special handling. For example, a monthly
average of a `NamedTuple` diagnostic accumulates and averages each component
independently, and then the averaged components are written as separate
variables, just like the instantaneous case.

## Why splitting happens at the writer

NetCDF (as used here, following the CF conventions) stores plain numeric
arrays with named dimensions; it has no notion of a "`NamedTuple` variable".
More importantly, the `NetCDFWriter` resamples every field onto a rectangular
grid using `ClimaCore`'s remapper, and the remapper operates on scalar fields
only.

For this reason a `NamedTuple` diagnostic is split into its scalar components
*before* remapping: each component is a zero-copy scalar view, remapped and
written through the ordinary scalar path. Everything upstream of the writer
(storage, accumulators, reductions, the averaging hook) already works on
`NamedTuple`-valued fields directly, so the split is confined to the writer
boundary and scalar diagnostics are completely unaffected.

The `HDF5Writer`/`DictWriter` do not remap, so they keep the field whole and
need no splitting.
