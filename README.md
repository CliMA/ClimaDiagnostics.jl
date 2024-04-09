<h1 align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="logo-white.svg">
  <source media="(prefers-color-scheme: light)" srcset="logo.svg">
  <img alt="Shows the logo of ClimaDiagnostics, with a globe and magnifying glasses" src="logo.svg" width="180px">
</picture>

ClimaDiagnostics.jl
</h1>

`ClimaDiagnostics.jl` provides a simple framework to add diagnostics to `CliMA`
simulations.

`ClimaDiagnostics.jl`defines two important concepts:
- `DiagnosticVariable`: A recipe to compute a diagnostic from the integrator
  alongside with names, units, comments.
- `ScheduledDiagnostic`: When to compute and output the `DiagnosticVariable` and
  what type of accumulation to perform.

To add the diagnostics to a simulation from a list of `ScheduledDiagnostic`s,
just redefine your integrator with `IntegratorWithDiagnostics(integrator,
scheduled_diagnostics)`.

The [documentation](https://clima.github.io/ClimaDiagnostics.jl/dev/) is rich
and comprehensive, please find all more information there.
