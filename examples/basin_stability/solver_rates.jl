using Graphs
using Plots
using OrdinaryDiffEq

using Revise
using MultilevelChainSampler

imgdir="$(@__DIR__)/imgs"

function fit_rate(X,Y)
    x = log.(X); y = log.(Y)
    A = hcat(x, ones(length(x)))
    Apinv = inv(A' * A)* A'
    z = Apinv * y
    a = z[1]
    a = round(a, digits=4)
    return a
end


N=24; K = 1.0; a = 0.5
G = complete_graph(N)
W = ones(nv(G)); W[1:nv(G)÷2] .= -1
g = PowerGrid(G,W,K,a)
timespan = (0.0, 20.0)

u_fix = synchronous_state(g)
perturbation = [ π/2, 10.0 ]
u₀ = u_fix[:]; u₀[ [1,N+1] ] += perturbation
problem = ODEProblem(swing_dynamics!, u₀, timespan, g)

# Plot high resolution dynamics
println("Solving with high resolution...")
timed = @elapsed fine_solution = solve(problem, Rodas5(), abstol=1e-12, reltol=1e-12);
print(); println("High resolution finished in $timed seconds.")
t_fine = fine_solution.t
u_fine = fine_solution.u[end]


function plot_rates(solvers, levels, values; labels="SOLVER, rate=RATE", ignore=0)

    plot(); xlims!(extrema(levels)...)
    for (i,solver) in enumerate(solvers)
        sol = string(solver)
        sol = split(string(sol),"(")[1]
        sol = split(string(sol),"{")[1]

        rate = fit_rate(levels[1+ignore:end], values[i,1+ignore:end])
        label = replace(labels, "SOLVER" => sol, "RATE" => "$rate")

        #plot!(levels, values[i,:]; label)
        plot!(levels[1+ignore:end], values[i,1+ignore:end]; label)
    end
    plot!(xscale=:log10, yscale=:log10, minorgrid=true)
    xlabel!("Level h")
    return Plots.current()
end


## Non-adaptive

solvers = [Euler(), ImplicitEuler(), RK4(), Tsit5()]
kwargs = (; adaptive=false)
timesteps = 0.1 .^ [ .1:.04:4 ... ]
errors = zeros(length(solvers), length(timesteps))
costs  = zeros(length(solvers), length(timesteps))
for (i,solver) in enumerate(solvers)
    for (j,h) in enumerate(timesteps)
        costs[i,j] = @elapsed solution = solve(problem, solver; dt=h, kwargs...)
        u = solution.u[end]
        errors[i,j] = abs(u[1] - u_fine[1]) / abs(u_fine[1])
    end
end

plot_rates(solvers, timesteps, errors; labels="SOLVER, α=RATE")
ylabel!("Relative error")
title!("Fixed-step Solvers - Convergence")
plot!(legend=:topleft)
plot!([extrema(timesteps)... ], repeat([1e-12],2), color="black", linestyle=:dash, label="")
savefig("$imgdir/solvers_fixed_error.png")
display(Plots.current())

plot_rates(solvers, timesteps, costs; labels="SOLVER, γ=RATE", ignore=3)
ylabel!("Evaluation time")
title!("Fixed-step Solvers - Costs")
savefig("$imgdir/solvers_fixed_costs.png")
display(Plots.current())


## Adaptive

solvers = [RK4(), Tsit5(), Rodas3(), Rodas5()]
kwargs = (; adaptive=true)
tolerances = 0.1 .^ [ 0:.05:7 ... ]
errors = zeros(length(solvers), length(tolerances))
costs  = zeros(length(solvers), length(tolerances))
for (i,solver) in enumerate(solvers)
    for (j,h) in enumerate(tolerances)
        costs[i,j] = @elapsed solution = solve(problem, solver; reltol=h, abstol=h,  kwargs...)
        u = solution.u[end]
        errors[i,j] = abs(u[1] - u_fine[1]) / abs(u_fine[1])
    end
end

plot_rates(solvers, tolerances, errors; labels="SOLVER, α=RATE")
ylabel!("Relative error")
title!("Adaptive solvers - Convergence")
plot!(legend=:topleft)
savefig("$imgdir/solvers_adaptive_error.png")
display(Plots.current())

plot_rates(solvers, tolerances, costs; labels="SOLVER, γ=RATE", ignore=3)
ylabel!("Evaluation time")
title!("Adaptive solvers - Costs")
savefig("$imgdir/solvers_adaptive_costs.png")
display(Plots.current())
