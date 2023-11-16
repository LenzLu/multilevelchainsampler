using Graphs

using Revise
using MultilevelChainSampler


clearconsole(); println(repeat("=",80)... )
println("Single node connected to clique ")

G = complete_graph(9)
add_vertex!(G); add_edge!(G, nv(G), 1)

W = ones(nv(G)); W[1:nv(G)÷2] .= -1
grid = PowerGrid(G,W,1.0,1.0)

t = @elapsed S = nodal_basin_stability(grid, nv(G))
println("Stability computation took $t seconds.")

println("Stability $S")
