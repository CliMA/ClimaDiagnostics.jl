using Test
import ClimaDiagnostics.DiagnosticVariables

@testset "DiagnosticVariable" begin
    # First, create a diagnostic variable

    # All the arguments
    var = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        standard_name = "my_test",
        units = "m",
        comments = "It works!",
        compute! = (out, u, p, t) -> 1,
    )

    @test DiagnosticVariables.short_name(var) == "my"

    # The minimum number of arguments required
    var =
        DiagnosticVariables.DiagnosticVariable(; compute! = (out, u, p, t) -> 1)

    @test DiagnosticVariables.short_name(var) == ""

end
