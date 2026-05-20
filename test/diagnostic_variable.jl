using Test
import Logging
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

@testset "component_units" begin
    cu = DiagnosticVariables.component_units

    # String: shared across every component (and the scalar/empty chain)
    @test cu("m/s", ()) == "m/s"
    @test cu("m/s", (:a,)) == "m/s"

    # NamedTuple: resolved per key (flat)
    @test cu((; a = "m", b = "kg"), (:a,)) == "m"
    @test cu((; a = "m", b = "kg"), (:b,)) == "kg"

    # NamedTuple: descended along a nested property chain
    @test cu((; warm = (; au = "1/s", ac = "1/s")), (:warm, :au)) == "1/s"

    # NamedTuple: flat-leaf fallback for a nested chain
    @test cu((; au = "1/s"), (:warm, :au)) == "1/s"

    # Function: called with the property chain
    @test cu(c -> join(string.(c), "."), (:warm, :au)) == "warm.au"

    # No entry returns `nothing`; an explicit "" is a deliberate empty value
    @test cu((; a = "m"), (:zzz,)) === nothing
    @test cu((; a = ""), (:a,)) == ""
    @test cu((; warm = (; au = "")), (:warm, :au)) == ""

    # An `Int` chain element (Tuple/NTuple-element index) must NOT positional-
    # index the units `NamedTuple` (regression: `haskey((; a=1), 1) == true`).
    @test cu((; a = "m", b = "kg"), (1,)) === nothing
    @test cu((; warm = (; au = "1/s")), (:warm, 1)) === nothing
    # ...but the flat-leaf fallback still resolves a `Symbol` leaf past an `Int`
    @test cu((; au = "1/s"), (2, :au)) == "1/s"
end

@testset "units_attribute" begin
    ua = DiagnosticVariables.units_attribute

    # String: returned unchanged, regardless of the property chains
    @test ua("kg/m3", [()]) == "kg/m3"
    @test ua("kg/m3", [(:a,), (:b,)]) == "kg/m3"

    # NamedTuple (flat): one "chain: units" entry per component
    @test ua((; cond = "x", evap = "y"), [(:cond,), (:evap,)]) ==
          "cond: x, evap: y"

    # NamedTuple (nested): chains are flattened with `.`
    @test ua(
        (; warm = (; au = "a", ac = "b"), tot = "c"),
        [(:warm, :au), (:warm, :ac), (:tot,)],
    ) == "warm.au: a, warm.ac: b, tot: c"

    # NamedTuple: flat-leaf fallback is honoured (unlike `string(nt)`)
    @test ua((; au = "1/s"), [(:warm, :au)]) == "warm.au: 1/s"

    # Function: resolved per component
    @test ua(chain -> "1/s", [(:warm, :au), (:tot,)]) ==
          "warm.au: 1/s, tot: 1/s"

    # A component with no matching entry gets an empty units value
    @test ua((; cond = "x"), [(:cond,), (:evap,)]) == "cond: x, evap: "

    # `Int` chain entries (Tuple-element indices) resolve to an empty value
    # rather than positional-indexing the units `NamedTuple`.
    @test ua((; a = "m"), [(1,), (:a,)]) == "1: , a: m"
end

@testset "units_warnings" begin
    uw = DiagnosticVariables.units_warnings

    # `String` is never checked (shared; the default `units` is `""`)
    @test_logs min_level = Logging.Warn uw("kg/m3", [(:a,), (:b,)])
    @test_logs min_level = Logging.Warn uw("", [(:a,), (:b,)])

    # `Function`: a non-`nothing` return (incl. "") is explicit -> no logs
    @test_logs min_level = Logging.Warn uw(c -> "1/s", [(:a,), (:b,)])
    @test_logs min_level = Logging.Warn uw(
        c -> c == (:b,) ? "" : "u",
        [(:a,), (:b,)],
    )
    # `Function` returning `nothing` for a component: one warning
    @test_logs (:warn,) uw(c -> c == (:b,) ? nothing : "u", [(:a,), (:b,)])

    # Fully specified NamedTuple: no logs
    @test_logs min_level = Logging.Warn uw((; a = "x", b = "y"), [(:a,), (:b,)])

    # An entry explicitly set to "" is deliberate: no logs
    @test_logs min_level = Logging.Warn uw((; a = ""), [(:a,)])
    @test_logs min_level = Logging.Warn uw(
        (; warm = (; au = "")),
        [(:warm, :au)],
    )

    # A component with no matching key: one warning (a units entry that does
    # not correspond to any component is NOT warned about)
    @test_logs (:warn,) uw((; a = "m"), [(:b,)])

    # Flat-leaf fallback resolves `:au`; only `:tot` has no entry: one warning
    @test_logs (:warn,) uw((; au = "1/s"), [(:warm, :au), (:tot,)])

    # `Int` chain entries (Tuple-element indices) cannot match a `NamedTuple`
    # entry and resolve to no units: one warning.
    @test_logs (:warn,) uw((; a = "m"), [(1,), (:a,)])
end

@testset "DiagnosticVariable non-string units" begin
    # `units` may be a NamedTuple (per-component) ...
    v_nt = DiagnosticVariables.DiagnosticVariable(;
        compute = (u, p, t) -> 1,
        short_name = "x",
        units = (; a = "m", b = "s"),
    )
    @test v_nt.units == (; a = "m", b = "s")

    # ... or a function of the property chain
    units_fn = chain -> "u"
    v_fn = DiagnosticVariables.DiagnosticVariable(;
        compute = (u, p, t) -> 1,
        short_name = "x",
        units = units_fn,
    )
    @test v_fn.units === units_fn
end
