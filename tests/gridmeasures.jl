using CairoMakie
using Statistics
using Graphs
using Random

using Revise
using MultilevelChainSampler

##
Random.seed!(0)
G = erdos_renyi(10, 20, seed=0); N = nv(G)
W = 20 .* rand(N); W .-= mean(W)
g = PowerGrid(G, W, 10.0, 0.01)

println("Synchronous state")
println("ϕ* = ", g.syncstate[1:N])
println("ω* = ", g.syncstate[N+1:2N])

println("r_uni = ", r_uni(g))
println("aspl = ", aspl(g))

println("Stability ", basin_stability(g))