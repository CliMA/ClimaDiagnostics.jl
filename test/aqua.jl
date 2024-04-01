using Test
using ClimaDiagnostics
using Aqua

@testset "Aqua tests" begin
    Aqua.test_undefined_exports(ClimaDiagnostics)
    Aqua.test_stale_deps(ClimaDiagnostics)
    Aqua.test_deps_compat(ClimaDiagnostics)
    Aqua.detect_ambiguities(ClimaDiagnostics; recursive = true)
    Aqua.test_piracies(ClimaDiagnostics)
end
