# Saving the diagnostics

Do not forget to close your writers to avoid file corruption!

## `NetCDFWriter`

The `NetCDFWriter` resamples the input `Field` to a rectangular grid and saves
the output to a NetCDF file.

The `NetCDFWriter` relies on the `Remappers` module in `ClimaCore` to
interpolate onto the rectangular grid. Horizontally, this interpolation is a
Lagrange interpolation, vertically, it is a linear. This interpolation is not
conservative. Also note that, the order of vertical interpolation drops to zero
in the first and last vertical elements of each column.

To create a `NetCDFWriter`, you need to specify the target `ClimaCore` `Space`
and the output directory where the files should be saved. By default, the
`NetCDFWriter` appends to existing files and create new ones if they do not
exist. The `NetCDFWriter` does not overwrite existing data and will error out if
existing data is inconsistent with the new one.

The output in the `NetCDFWriter` roughly follows the CF conventions.

Each `ScheduledDiagnostic` is output to a different file with name determined by
calling the `output_short_name` on the `ScheduledDiagnostic`. Typically, these
files have names like `ta_1d_max.nc`, `ha_20s_inst.nc`, et cetera. The files
define their dimensions (`lon`, `lat`, `z`, ...). Time is always the first
dimension is any dataset.

Variables are saved as datasets with attributes, where the attributes include
`long_name`, `standard_name`, `units`...

```@docs
ClimaDiagnostics.Writers.NetCDFWriter
ClimaDiagnostics.Writers.write_field!
Base.close
```

## `DictWriter`

The `DictWriter` is a in-memory writer that is particularly useful for
interactive work and debugging.

```@docs
ClimaDiagnostics.Writers.DictWriter
ClimaDiagnostics.Writers.write_field!
```

## `HDF5Writer`

 The `HDF5Writer` writes the `Field` directly to an HDF5 file in such a way that
it can be later read and imported using the `InputOutput` module in `ClimaCore`.

The `HDF5Writer` writes one file per variable per timestep. The name of the file
is determined by the `output_short_name` field of the `ScheduledDiagnostic` that
is being output.

> Note: The `HDF5Writer` in `ClimaDiagnostics` is currently the least developed
> one. If you need this writer, we can expand it.

```@docs
ClimaDiagnostics.Writers.HDF5Writer
ClimaDiagnostics.Writers.write_field!
```
