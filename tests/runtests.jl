using Test
using MultilevelChainSampler

clearconsole(); println(repeat("=",80)... )
println("Testing... ")

@testset "Minimal" begin include("mhchain_gauss.jl") end
