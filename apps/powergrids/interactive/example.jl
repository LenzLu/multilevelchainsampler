include("nodal.jl")

using Graphs 
using Statistics
using Random 

##
Random.seed!(42)
G = erdos_renyi(10, 20, seed=0); N = nv(G)
W = 20 .* rand(N); W .-= mean(W)
g = PowerGrid(G, W, 10.0, 0.1)

println("Synchronous state")
println("ϕ* = ", g.syncstate[1:N])
println("ω* = ", g.syncstate[N+1:2N])

## 
fig = analyze_nodal_stability(g, 9, 128; tolerance=0.1, t_final=128.0)
display(fig)