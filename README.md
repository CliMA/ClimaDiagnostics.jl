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

`ClimaDiagnostics.jl`defined two important concepts:
- `DiagnosticVariable`: A recipe to compute a diagnostic from the integrator
  alongside with names, units, comments.
- `ScheduledDiagnostic`: When to compute and output the `DiagnosticVariable` and
  what type of accumulation to perform.

To add the diagnostics to a simulation from a list of `ScheduledDiagnostic`s,
one first needs to initialize them with `diagnostic_handler =
DiagnosticHandler(diangostics, Y, p, t; dt)`, and finally add the
`DiagnosticsCallback(diagnostic_handler)` callback to the integrator.

> :warning: README under construction. While we work on the README, you can find
> all the relevant information in the
> [documentation](https://clima.github.io/ClimaDiagnostics.jl/dev/).
