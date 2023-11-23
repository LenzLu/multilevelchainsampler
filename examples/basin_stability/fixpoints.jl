using Graphs
using GraphRecipes

using Revise
using MultilevelChainSampler


N=24; K = 1.0; a = 0.5
G = path_graph(N)
W = ones(nv(G)); W[1:nv(G)÷2] .= -1
g = PowerGrid(G,W,K,a)

u_fix = synchronous_state(g)
is_stable(g, u_fix)
