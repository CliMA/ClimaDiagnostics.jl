# Public APIs

## `ClimaDiagnostics`

```@docs
ClimaDiagnostics.IntegratorWithDiagnostics
```

## `Schedules`

```@docs
ClimaDiagnostics.Schedules.AbstractSchedule
ClimaDiagnostics.Schedules.short_name
ClimaDiagnostics.Schedules.long_name
ClimaDiagnostics.Schedules.DivisorSchedule
ClimaDiagnostics.Schedules.EveryStepSchedule
ClimaDiagnostics.Schedules.EveryDtSchedule
ClimaDiagnostics.Schedules.EveryCalendarDtSchedule
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
ClimaDiagnostics.Writers.interpolate_field!
ClimaDiagnostics.Writers.write_field!
Base.close
```
