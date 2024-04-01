using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "Format" begin @time include("format.jl") end
#! format: on

return nothing
