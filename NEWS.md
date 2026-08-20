# NEWS

v0.3.7
------

## Features

### `CalendarIntervalSchedule`: calendar-partitioned output windows

A new output schedule, `CalendarIntervalSchedule`, partitions time into exact calendar
windows of a given `Dates.Period` (e.g. `Dates.Month(1)`) and closes the window
`[date_last, date_last + period)` when the date reaches a period boundary. Unlike
`EveryCalendarDtSchedule`, it advances `date_last` to the exact calendar boundary (not the
overshooting step time) and reports the closed window to the writer.

Windows are exact left-closed partitions: a sample at a boundary belongs to the window
*starting* there and is folded after the previous window's output. As a consequence,
restarting a simulation exactly at a closed boundary (with `date_last` set to it) reproduces
an uninterrupted run exactly.

```julia
import ClimaDiagnostics.Schedules: CalendarIntervalSchedule
import Dates
sched = CalendarIntervalSchedule(Dates.Month(1); start_date = Dates.DateTime(2024, 1, 1))
```

An optional `spinup` date marks windows that *start* before it as inactive: no samples are
accumulated and nothing is written for them. This is convenient for excluding a model spin-up
period from a monthly/seasonal mean.

```julia
sched = CalendarIntervalSchedule(
    Dates.Month(1);
    start_date = Dates.DateTime(2024, 1, 1),
    spinup = Dates.DateTime(2024, 2, 1),   # [Jan 1, Feb 1) is skipped
)
```

### `IntervalSchedule`: wrap any schedule into interval windows

`IntervalSchedule` wraps an existing time-based schedule (e.g. `EveryDtSchedule`) and turns
its firing times into window boundaries: if it fires at `t1, t2, ...`, the windows are
`[t0, t1)`, `[t1, t2)`, ... with the same exact sample attribution, `spinup` gating, and
exact `time_bnds`/`date_bnds` as `CalendarIntervalSchedule`. Boundaries are the *realized*
firing times (they follow the stepper under overshoot); use `CalendarIntervalSchedule` when
the edges must lie on a calendar grid. The wrapped schedule is consulted by the wrapper and
must not be shared or called elsewhere, and it must only read `integrator.t` (so
step-counting schedules like `EveryStepSchedule` cannot be wrapped).

```julia
import ClimaDiagnostics.Schedules: IntervalSchedule, EveryDtSchedule
sched = IntervalSchedule(
    EveryDtSchedule(30 * 86400.0);
    start_date = Dates.DateTime(2024, 1, 1),
)
```

### New schedule protocol: `should_accumulate` / `output_period_bounds`

Schedules may now optionally implement two methods (both have generic fallbacks, so existing
schedules are unchanged):

- `should_accumulate(schedule, integrator)` — whether the current time is inside the active
  accumulation window currently open (default `true`). It is read before and after the
  schedule's output call; a sample that only belongs to the newly opened window (a boundary
  sample) is folded after the output.
- `output_period_bounds(schedule)` — the `[lo, hi)` window the schedule last closed, or
  `nothing` (default `nothing`).

## Behavior change (opt-in)

For time-reduction diagnostics whose output schedule reports a window via
`output_period_bounds` (currently only `CalendarIntervalSchedule`), the NetCDF writer now
stamps the schedule's **exact** window in `time_bnds`/`date_bnds`, timestamped at the window
start, instead of reconstructing `[previous_output, current]`. Diagnostics using any other
schedule are unaffected.

v0.3.6
------

- `IntegratorWithDiagnostics` no longer reads `integrator.callback.continuous_callbacks`.
  `CallbackSet` is now constructed from discrete callbacks only, matching the
  corresponding upcoming change in ClimaTimeSteppers v1.0.

v0.3.5
------

#### Allow ClimaTimeSteppers v0.9 PR[#168](https://github.com/CliMA/ClimaDiagnostics.jl/pull/168)

v0.3.4
------

#### Use updated ClimaTimeSteppers interface PR[#164](https://github.com/CliMA/ClimaDiagnostics.jl/pull/164)

#### Use Julia v1.12 in CI PR[#158](https://github.com/CliMA/ClimaDiagnostics.jl/pull/158)

#### Bump GHA versions PR[#163](https://github.com/CliMA/ClimaDiagnostics.jl/pull/163) PR[#162](https://github.com/CliMA/ClimaDiagnostics.jl/pull/162)

v0.3.3
------

## Features

### Support bilinear horizontal remapping

Adds support for a `horizontal_method` kwarg in the NetCDFWriter, to permit
bilinear remapping using ClimaCore's AbstractRemappingMethod.

v0.3.1
-------

## Features

### Support 1D column diagnostics

Added support for outputting the 1-dimensional FiniteDifferenceSpace diagnostics
when writer is initialized from 3D space.

### Add global attributes to NetCDF files

With this release, you can now add global attributes that are the same across
all NetCDF files for a given `NetCDFWriter`. For example, you may be interested
in specifying the `source` and `experiment` which are the same across all
NetCDF files produced for a single simulation. You can do this via the
`global_attribs` keyword argument for the `NetCDFWriter`. The global attributes
must be a subtype of an `AbstractDict{String, String}`. If the order of the
attributes matter, you may want to use an `OrderedDict` from
`OrderedCollections`.

```julia
writer = NetCDFWriter(
    space,
    output_dir;
    global_attribs = Dict("source" => "CliMA Coupler Simulation", "experiment" => "AMIP"),
)
```

## Output diagnostics in pressure coordinates

You can now output diagnostics in pressure coordinates. To do this, create an
instance of a `RealPressureLevelsMethod` with the pressure field and the current
simulation time. Then, pass the `RealPressureLevelsMethod` to the
`z_sampling_method` keyword argument to the `NetCDFWriter`.

```julia
z_sampling_method = ClimaDiagnostics.Writers.RealPressureLevelsMethod(
            pressure_field,
            t,
        )
netcdf_writer = CAD.NetCDFWriter(
        ClimaDiagnostics.Writers.pressure_space(z_sampling_method),
        output_dir,
        num_points = (360, 180, 10); # the number of vertical points (10) is ignored
        sync_schedule = CAD.EveryStepSchedule(),
        z_sampling_method,
    )
```

v0.3.0
-------

## Breaking changes

- ![][badge-💥breaking]  Reduced NetCDF diagnostics are now timestamped at the start
  of the reduction period instead of the end. Instantaneous diagnostics are unchanged.

## Specify points for horizontal interpolation

With this relase, `horizontal_pts` is added as a keyword argument to
`NetCDFWriter` which allow users to specify specific points for interpolation in
the horizontal space.

In the example below, interpolation is done at the points along the longitudes
from 0.5 degrees to 179.5 degrees and the latitudes from -90.0 degrees to 90.0
degrees with a step of 1 degree for both.

```julia
writer = NetCDFWriter(
    space, # 2D space with longitudes and latitudes
    output_dir;
    horizontal_pts = (0.5:179.5, -90.0:90.0),
)
```

## Bug fixes

Fixed `num_points` representing Lat-Long-Z for a box domain.

v0.2.14
-------

## Bug fixes

Fixed the default target coordinates for the `NetCDFWriter` for boxes with only
one point in the horizontal space.

v0.2.13
-------

## Features

### Add time and date bounds for netCDF files

In a netCDF file produced by ClimaDiagnostics, there are now the dimensions
`time_bnds` and `date_bnds`. Each time or date value is a representative of the
corresponding time or date bound. For example, if the time is `10.0` and
corresponding time bound is `[0.0, 10.0]`, then the time of `10.0` represents
the interval `[0.0, 10.0]`. If one knows that the data represents a time
average, then the time of `10.0` is the time average over the interval
`[0.0, 10.0]`.

### NetCDF writer now defaults to a reasonable number of points

`ClimaDiagnostics.Writers.NetCDF` now has a new default argument that depends on
the input Space. With this new default, obtained by calling the
`ClimaDiagnostics.Writers.default_num_points(space)` function, the output
diagnostics will be sampled with approximately the same resolution as the given
`space`.

### Support for `lazy`

Starting version `0.2.13`, `ClimaDiagnostics` supports diagnostic variables
specified with un-evaluated expressions (as provided by
[LazyBroadcast.jl](https://github.com/CliMA/LazyBroadcast.jl)).

Instead of

```julia
function compute_ta!(out, state, cache, time)
    if isnothing(out)
        return state.ta .- 273.15
    else
        out .= state.ta .- 273.15
    end
end
```

You can now write

```julia
import LazyBroadcast: lazy

function compute_ta(state, cache, time)
    return lazy.(state.ta .- 273.15)
end
```

Or, for `Field`s

```julia
function compute_ta(state, cache, time)
    return state.ta
end
```

v0.2.12
-------

## Features

### Support with ITime

This release adds support for working with `ITime`s. In particular,
`EveryCalendarDtSchedule` is compatible with `ITime` and a new constructor is
provided to make an `EveryCalendarDtSchedule` using `ITime`s. Lastly, there are
small changes to the writers to support `ITime`s.

v0.2.12
-------

## Bug fixes

- `NetCDFWriter` now correctly writes purely vertical and point spaces.

v0.2.11
-------

## Bug fixes

- Times in `DictWriter` are now correctly sorted.

v0.2.10
-------

## Bug fixes

Fixed broken `start_date` feature.

v0.2.9
-------

## Features

### Add a `start_date` attribute to NetCDFWriter PR [#94](https://github.com/CliMA/ClimaDiagnostics.jl/pull/94)

Prior to this version, users had to go to the simulation to find the start date.
It can now be saved as an attribute, making it easily accessible.
To do so, users need to pass the kwarg `start_date` when calling `NetCDFWriter`.

## Bug fixes

### Acquiring ownership with `compute!` PR [#88](https://github.com/CliMA/ClimaDiagnostics.jl/pull/88)

Prior to this version, `ClimaDiagnostics` would directly store use the output
returned by `compute!` functions the first time they are called. This leads to
problems when the output is a reference to an existing object since multiple
diagnostics would modify the same object. Now, `ClimaDiagnostics` makes a copy
of the return object so that it is no longer necessary to do so in the
`compute!` function.

### Correctly de-duplicate `ScheduledDiagnostics` [#93](https://github.com/CliMA/ClimaDiagnostics.jl/pull/93)

This version fixes a bug where `ScheduledDiagnostics` were not correctly
de-duplicated because `==` was not implemented correctly.

v0.2.8
-------

## Bug fixes

- `IntegratorWithDiagnostics` advertised a feature that was not implemented:
  `IntegratorWithDiagnostics` claimed that passing `state_name` and `cache_name`
  would allow users to customize the name of the state and cache inside the
  integrator. Now, this is implemented.

v0.2.7
-------

## Bug fixes

- `scheduled_diagnostics` are now internally saved as vectors instead of tuples.
  This has significant compile-time/inference benefits.

v0.2.6
-------

## Features

### More matadata in NetCDF files

Release `0.2.6` improves compatibility with CF conventions by adding

- standard and long name for the `time`, `longitude`, and `latitude` dimensions

## Bug fixes

- The default constructor for `ScheduleDiagnostic`s no longer uses reference of
  `Schedule`s but create a new copy.

### Deprecations

`reference_date` was renamed to `start_date` and `t_start` was dropped from the
constructors for the schedules. These changes are due to the fact that these
arguments should not be needed.

v0.2.5
-------

## Features

### Add support for box spaces with LatLong points in `NetCDFWriter`

The `NetCDFWriter` can now work with regional boxes with `LatLong` points. Due
to incompatibility in `ClimaCore`, only `LatLong` points are supported (and not
`LongLat` points). This means that the box has to be created with latitude on
the x axis and longitude on the y axis.

## Bug fixes

- Ensure that `DictWriter` can only be constructed with dictionary-like objects.

v0.2.4
-------

- Add `EveryCalendarDtSchedule` for schedules with calendar periods.

v0.2.3
-------

- Detect and ignore duplicated diagnostics.

v0.2.2
-------

- Fix support for caches that are not NamedTuples in `NetCDFWriter`.

v0.2.1
-------

- Fix support for purely horizontal spaces in `NetCDFWriter`.

v0.2.0
-------

- ![][badge-💥breaking] The `NetCDFWriter` now outputs points on the vertical levels by default.
- ![][badge-💥breaking] `disable_vertical_interpolation` is removed in favor of `z_sampling_method`.

[badge-💥breaking]: https://img.shields.io/badge/💥BREAKING-red.svg
