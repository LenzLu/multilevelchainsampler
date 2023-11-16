using Graphs
using Plots
using OrdinaryDiffEq

using Revise
using MultilevelChainSampler



#grid = PowerGrid( path_graph(2), [1,-1], 1.0, 1.0)
G = path_graph(10); add_edge!(G, 10, 1)
W = ones(nv(G)); W[nv(G)÷2:end] .= -1
grid = PowerGrid(G,W,1.0,1.0)

u_fix = steady_state(grid)
#reshape([(2*rand()-1)*π, randn()],1,2),
time_steps = 10 .^ [ -4.0:0.5:0.0 ... ]
algos = [ :RK4, :Tsit5, :Heun ]
times = Dict()
stabs = Dict()
for algo in algos
    solver = eval(algo)()
    times[algo] = zeros(length(time_steps))
    stabs[algo] = zeros(length(time_steps))
    for (i,dt) in enumerate(time_steps)
        T = @elapsed s = basin_stability( grid; n_sample=100, u_fix, solver, adaptive=false, dt )
        times[algo][i] = T
        stabs[algo][i] = s
    end
end


plot(); xlims!(extrema(time_steps)...)
for algo in algos
    name = string(algo)
    plot!(time_steps, times[algo], label="Solver $name ")
end
xlabel!("time step")
ylabel!("computation time")
xaxis!(:log10); yaxis!(:log10)
display(Plots.current())

plot(); xlims!(extrema(time_steps)...)
for algo in algos
    name = string(algo)
    plot!(time_steps, stabs[algo], label="Solver $name ")
end
xlabel!("time step")
ylabel!("stability")
xaxis!(:log10)
display(Plots.current())
