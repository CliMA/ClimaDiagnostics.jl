# Saving the diagnostics


Do not forget to close your writers to avoid file corruption!

## `NetCDFWriter`

```@docs
ClimaDiagnostics.Writers.NetCDF
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

The `HDF5Writer` in `ClimaDiagnostics` is currently the least developed one.

```@docs
ClimaDiagnostics.Writers.NetCDF
ClimaDiagnostics.Writers.write_field!
Base.close
```
