
using Graphs 
using Statistics
using Random 

using Revise
using MultilevelChainSampler

##
Random.seed!(0)
G = erdos_renyi(10, 20, seed=0); N = nv(G)
W = 20 .* rand(N); W .-= mean(W)
g = PowerGrid(G, W, 10.0, 0.1)

println("Synchronous state")
println("ϕ* = ", g.syncstate[1:N])
println("ω* = ", g.syncstate[N+1:2N])

## 
fig = analyse_nodal_stability(g, 3, 256; tolerance=0.1, t_final=128.0)
