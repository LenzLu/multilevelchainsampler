using Graphs
using Plots

using Revise
using MultilevelChainSampler

clearconsole(); println(repeat("=",80)... )
println("Estimating convergence rates ")
imgdir = "$(@__DIR__)/imgs"

N=8; K = 5.0; a = 0.5
G = complete_graph(N)
W = ones(nv(G)); W[1:nv(G)÷2] .= -1
grid = PowerGrid(G,W,K,a)
timespan = (0.0, 20.0)

u_fix = synchronous_state(grid)
perturbation = [ π/2, 10.0 ]
u₀ = u_fix[:]; u₀[ [1,N+1] ] += perturbation
problem = ODEProblem(swing_dynamics!, u₀, timespan, grid)

# Plot high resolution dynamics
println("Solving with high resolution...")
timed = @elapsed fine_solution = solve(problem, Rodas5(), abstol=1e-19, reltol=1e-19)
println("Solver finished in $timed seconds.");
t_fine = fine_solution.t
u_fine = hcat(fine_solution.u ... )

colors = [ RGB(0, i/N, 1-i/N) for i=1:N ]
plot(); xlims!(timespan...)
for i = 1:N
    plot!(t_fine, u_fine[i,:], label="\$ϕ_$i\$", c=colors[i])
end
xlabel!("time"); ylabel!("angle")
savefig("$imgdir/trajectory_angle.png")
display(Plots.current())

plot(); xlims!(timespan...)
for i=1:N
    plot!(t_fine, u_fine[i+N,:], label="\$ω_$i\$", c=colors[i])
end
xlabel!("time"); ylabel!("frequency")
savefig("$imgdir/trajectory_frequency.png")
display(Plots.current())

plot(); xlims!(timespan...)
plot!(fine_solution.t, u_fine[1,:], label="fine solution")
for h = .5 .^[2:5 ... ]
    solution = solve(problem, RK4(), adaptive=false, dt=h)
    u = hcat(solution.u ... )
    plot!(solution.t, u[1,:], label="RK4 Δt=$h")
end
xlabel!("time"); ylabel!("angle")
savefig("$imgdir/trajectories_angle.png")
display(Plots.current())
