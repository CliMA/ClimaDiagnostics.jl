using Documenter
import ClimaDiagnostics

@testset "Test docstrings" begin
    doctest(ClimaDiagnostics; manual = false)
end
