# Saving the diagnostics

Writers are needed to save the computed diagnostics.

`ClimaDiagnostics` comes with three writers:
- `NetCDFWriter`, to interpolate and save to NetCDF files;
- `DictWriter`, to save `Field`s to dictionaries in memory;
- `HDF5Writer`, to save `Field`s to HDF5 files.

(There is an additional `DummyWriter` that does nothing. It is mostly used
internally for testing and debugging.)

Users are welcome to implement their own writers. A writer has to be a subtype
of `AbstractWriter`and has to implement the `interpolate_field!` and
`write_field!` methods. `interpolate_field!` can return `nothing` is no
interpolation is needed.

## `NetCDFWriter`

The `NetCDFWriter` resamples the input `Field` to a rectangular grid and saves
the output to a NetCDF file.

The `NetCDFWriter` relies on the `Remappers` module in `ClimaCore` to
interpolate onto the rectangular grid. Horizontally, this interpolation is a
Lagrange interpolation, vertically, it is a linear. This interpolation is not
conservative. Also note that, the order of vertical interpolation drops to zero
in the first and last vertical elements of each column.

To create a `NetCDFWriter`, you need to specify the source `ClimaCore` `Space`
and the output directory where the files should be saved. By default, the
`NetCDFWriter` appends to existing files and create new ones if they do not
exist. The `NetCDFWriter` does not overwrite existing data and will error out if
existing data is inconsistent with the new one.

Optionally (recommended), you can pass an optional argument `start_date`, which
will be saved as an attribute of your NetCDF file, easily accessible.

`NetCDFWriter`s take as one of the inputs the desired number of points along
each of the dimensions. For the horizontal dimensions, points are sampled
linearly. For the vertical dimension, the behavior can be customized by passing
the `z_sampling_method` variable. When `z_sampling_method =
ClimaDiagnostics.Writers.LevelMethod()`, points evaluated on the grid levels
(and the provided number of points ignored), when `z_sampling_method =
ClimaDiagnostics.Writers.FakePressureLevelsMethod()`, points are sampled
uniformly in simplified hydrostatic atmospheric model.

The output in the `NetCDFWriter` roughly follows the CF conventions.

Each `ScheduledDiagnostic` is output to a different file with name determined by
calling the `output_short_name` on the `ScheduledDiagnostic`. Typically, these
files have names like `ta_1d_max.nc`, `ha_20s_inst.nc`, et cetera. The files
define their dimensions (`lon`, `lat`, `z`, ...). Time is always the first
dimension is any dataset.

Do not forget to close your writers to avoid file corruption!

To help reducing data loss, `NetCDFWriter` can force __syncing__, i.e. flushing
the values to disk. Usually, NetCDF buffers writes to disk (because they are
expensive), meaning values are not immediately written but are saved to disk in
batch. This can result in data loss, and it is often useful to force NetCDF to
write to disk (this is especially the case when working with GPUs). To do so,
you can pass the `sync_schedule` function to the constructor of `NetCDFWriter`.
When not `nothing`, `sync_schedule` is a callable that takes one argument (the
`integrator`) and returns a bool. When the bool is true, the files that were
modified since the last `sync` will be `sync`ed. For example, to force sync
every 1000 steps, you can pass the
`ClimaDiagnostics.Schedules.DivisorSchedule(1000)` schedule. By default, on
GPUs, we call `sync` at the end of every time step for those files that need to
be synced.

Variables are saved as datasets with attributes, where the attributes include
`long_name`, `standard_name`, `units`...

Global attributes can be added to the NetCDF files via the `global_attribs`
keyword argument for the `NetCDFWriter`. For example, you may want to specify
the `source` and `experiment` attributes, which are the same across all NetCDF
files produced for a single simulation.

```julia
writer = NetCDFWriter(
    space, # 2D space with longitudes and latitudes
    output_dir;
    global_attribs = Dict("source" => "CliMA Coupler Simulation", "experiment" => "AMIP"),
)
```

The global attributes must be a subtype of `AbstractDict{String, String}`. If
the order of the attributes matters, you may want to use an
[`OrderedDict`](https://juliacollections.github.io/OrderedCollections.jl/dev/#OrderedDicts)
from `OrderedCollections.jl`.

!!! note
    The `NetCDFWriter` cannot save raw `ClimaCore.Fields`, only fields that are
    resampled onto a Cartesian grids are supported. If you need such capability,
    consider using the [`ClimaDiagnostics.Writers.HDF5Writer`](@ref).

```@docs; canonical = false
ClimaDiagnostics.Writers.NetCDFWriter(space::Spaces.AbstractSpace, output_dir)
ClimaDiagnostics.Writers.interpolate_field!(writer::ClimaDiagnostics.Writers.NetCDFWriter, field, diagnostic, u, p, t)
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.NetCDFWriter, field, diagnostic, u, p, t)
ClimaDiagnostics.Writers.sync(writer::ClimaDiagnostics.Writers.NetCDFWriter)
Base.close(writer::ClimaDiagnostics.Writers.NetCDFWriter)
```

Sampling methods for the vertical direction:
```@docs
ClimaDiagnostics.Writers.AbstractZSamplingMethod
ClimaDiagnostics.Writers.LevelsMethod
ClimaDiagnostics.Writers.FakePressureLevelsMethod
```

### Output diagnostics in pressure coordinates

To write diagnostics in pressure coordinates, you must pass a
`RealPressureLevelsMethod` to the `NetCDFWriter` which specifies that
diagnostics should be interpolated to pressure levels for the vertical
direction.

To create a `RealPressureLevelsMethod`, you must pass a pressure field and the
current simulation time. The `RealPressureLevelsMethod` instance is then passed
to `NetCDFWriter` via the `z_sampling_method` keyword argument, which causes all
diagnostics to be output in pressure coordinates.

```julia
z_sampling_method = ClimaDiagnostics.Writers.RealPressureLevelsMethod(
            pressure_field,
            t,
            pressure_attribs = (; units = "Pa"),
            pressure_levels = [0.0, 10000.0]
        )
netcdf_writer = NetCDFWriter(
        ClimaDiagnostics.Writers.pressure_space(z_sampling_method),
        output_dir,
        num_points = (360, 180, 10); # the number of vertical points (10) is ignored
        sync_schedule = CAD.EveryStepSchedule(),
        z_sampling_method,
    )
```

In the example above, a `RealPressureLevelsMethod` instance is constructed. The
keyword argument `pressure_attribs` is a `NamedTuple` that populates the
attributes of the `pressure_level` dimension in the NetCDF file. The keyword
argument `pressure_levels` is a vector of sorted pressure levels that are being
interpolated to. When using `RealPressureLevelsMethod`, the space passed to the
`NetCDFWriter` must be a space with pressure as the vertical coordinate. You can
obtain this space by calling `pressure_space` on a `RealPressureLevelsMethod`
instance.

When using `RealPressureLevelsMethod`, the `NetCDFWriter` appends `_pressure` to
the output filename.

For more detail about how pressure interpolation is done, see the ClimaCore
[documentation](https://clima.github.io/ClimaCore.jl/dev/remapping/#Interpolating-to-pressure-coordinates).

```@docs; canonical = false
ClimaDiagnostics.Writers.RealPressureLevelsMethod
ClimaDiagnostics.Writers.RealPressureLevelsMethod(pfull_field, t)
```

## `DictWriter`

The `DictWriter` is a in-memory writer that is particularly useful for
interactive work and debugging.

```@docs
ClimaDiagnostics.Writers.DictWriter()
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.DictWriter, field, diagnostic, u, p, t)
```

## `HDF5Writer`

 The `HDF5Writer` writes the `Field` directly to an HDF5 file in such a way that
it can be later read and imported using the `InputOutput` module in `ClimaCore`.

The `HDF5Writer` writes one file per variable per timestep. The name of the file
is determined by the `output_short_name` field of the `ScheduledDiagnostic` that
is being output.

> Note: The `HDF5Writer` in `ClimaDiagnostics` is currently the least developed
> one. If you need this writer, we can expand it.

```@docs; canonical=false
ClimaDiagnostics.Writers.HDF5Writer
ClimaDiagnostics.Writers.write_field!(writer::ClimaDiagnostics.Writers.HDF5Writer, field, diagnostic, u, p, t)
Base.close(writer::ClimaDiagnostics.Writers.HDF5Writer)
```
