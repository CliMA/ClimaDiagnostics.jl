using Documenter
using ClimaDiagnostics

pages = ["Overview" => "index.md"]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)
format = Documenter.HTML(
    prettyurls = !isempty(get(ENV, "CI", "")),
    collapselevel = 1,
    mathengine = mathengine,
)

makedocs(
    sitename = "ClimaDiagnostics.jl",
    authors = "Gabriele Bozzola",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    strict = false,
    clean = true,
    modules = [ClimaDiagnostics],
)

deploydocs(
    repo = "github.com/CliMA/ClimaDiagnostics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
