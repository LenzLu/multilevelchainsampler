
using Graphs 
using Statistics

using Revise
using MultilevelChainSampler


G = erdos_renyi(50, 120) 
W = ones(nv(G)); W[1+nv(G)÷2:end] .= -1
g = PowerGrid(G, W, 6.0, 0.01)

analyse_nodal_stability(g, 1)