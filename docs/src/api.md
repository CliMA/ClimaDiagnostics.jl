# Public APIs

## `ClimaDiagnostics`

```@docs
ClimaDiagnostics.DiagnosticsHandler
ClimaDiagnostics.DiagnosticsCallback
```

## `Callbacks`

```@docs
ClimaDiagnostics.Callbacks.Callback
ClimaDiagnostics.Callbacks.CallbackOrchestrator
ClimaDiagnostics.Callbacks.AbstractSchedule
ClimaDiagnostics.Callbacks.short_name
ClimaDiagnostics.Callbacks.long_name
ClimaDiagnostics.Callbacks.DivisorSchedule
ClimaDiagnostics.Callbacks.EveryStepSchedule
ClimaDiagnostics.Callbacks.EveryDtSchedule

```

## `DiagnosticVariables`

```@docs
ClimaDiagnostics.DiagnosticVariables.DiagnosticVariable
ClimaDiagnostics.DiagnosticVariables.short_name
ClimaDiagnostics.DiagnosticVariables.long_name
ClimaDiagnostics.DiagnosticVariables.descriptive_short_name
ClimaDiagnostics.DiagnosticVariables.descriptive_long_name
```

## `ScheduledDiagnostics`

```@docs
ClimaDiagnostics.ScheduledDiagnostics.ScheduledDiagnostic
ClimaDiagnostics.ScheduledDiagnostics.output_short_name
ClimaDiagnostics.ScheduledDiagnostics.output_long_name
```


## `Writers`

```@docs
ClimaDiagnostics.AbstractWriter
ClimaDiagnostics.Writers.DictWriter
ClimaDiagnostics.Writers.NetCDFWriter
ClimaDiagnostics.Writers.HDF5Writer
ClimaDiagnostics.Writers.write_field!
Base.close
```
