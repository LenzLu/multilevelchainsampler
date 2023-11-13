include("../../PowergridStability/all.jl")

using Statistics

K=9.0; α=1.0
G = path_graph(10)
W = rand(nv(G)); W .-= mean(W)
grid = PowerGrid(G, W, K, α)

u = steady_state(grid)
solve_system(grid, u)
