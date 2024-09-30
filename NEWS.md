# NEWS

main
-------

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

### Add support for box spaces with LatLong points in `NetCDFWriter`.

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
