project = dirname(@__DIR__) * "/src"
include("$project/ChainSamplers/sampler.jl")
include("$project/NetworkEnsembles/all.jl")
include("$project/PowerGrids/module.jl")


## 
G = erdos_renyi(10, 20, seed=0); N = nv(G)
g = PowerGrid(G)

println("Synchronous state")
println("ϕ* = ", g.syncstate[1:N])
println("ω* = ", g.syncstate[N+1:2N])

plt = analyse_network(g, true)


