using Graphs

using Revise
using MultilevelChainSampler


clearconsole(); println(repeat("=",80)... )
println("Single node connected to clique ")

K = 9.0; a = 0.5
for N = 2:2:10
    println("N=$N nodes")

    G = complete_graph(N-1)
    add_vertex!(G); add_edge!(G, N, 1)

    W = ones(nv(G)); W[1:nv(G)÷2] .= -1
    grid = PowerGrid(G,W,K,a)

    @time S,_ = nodal_basin_stability(grid, 1 )
    #println("Stability computation took $t seconds.")

    println("Stability $S")
end
