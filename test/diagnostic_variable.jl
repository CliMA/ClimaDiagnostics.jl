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

    # Passing both compute and compute!
    @test_throws ErrorException DiagnosticVariables.DiagnosticVariable(;
        compute! = (out, u, p, t) -> 1,
        compute = (u, p, t) -> 1,
    )

    # Passing no compute or compute!
    @test_throws ErrorException DiagnosticVariables.DiagnosticVariable(;)
end

@testset "DiagnosticVariable show" begin
    var = DiagnosticVariables.DiagnosticVariable(;
        short_name = "my",
        long_name = "My test",
        compute! = (out, u, p, t) -> 1,
    )

    out = sprint(show, MIME("text/plain"), var)
    @test occursin("DiagnosticVariable", out)
    @test count(==('\n'), out) <= 10

    out2 = sprint(show, var)
    @test occursin("DiagnosticVariable", out2)
    @test !occursin('\n', out2)

    out3 = sprint(show, MIME("text/plain"), var; context = :compact => true)
    @test out2 == out3

    out_summary = sprint(summary, var)
    @test occursin("DiagnosticVariable", out_summary)
    @test !occursin('\n', out_summary)
end
