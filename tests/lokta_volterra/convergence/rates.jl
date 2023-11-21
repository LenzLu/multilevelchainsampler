include("../dynamics.jl")

using Plots
using Revise
using MultilevelChainSampler

imgdir="$(@__DIR__)/../imgs"




problem = ODEProblem(lokta_volterra!, u₀, (0,t_final), p)

print("Solving with high resolution...")
timed = @elapsed fine_solution = solve(problem, Rodas5(); adaptive=true, abstol=1e-16, reltol=1e-16)
print(); println("High resolution finished in $timed seconds.")
t_fine = fine_solution.t
u_fine = fine_solution.u[end]


function fit_rate(X,Y)
    x = log.(X); y = log.(Y)
    A = hcat(x, ones(length(x)))
    Apinv = inv(A' * A)* A'
    z = Apinv * y
    a = z[1]
    a = round(a, digits=4)
    return a
end

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


solvers = [Euler(), ImplicitEuler(), Heun(), RK4(), Tsit5(), Rodas5()]
kwargs = (; adaptive=false)
timesteps = 0.5 .^ [ 0:0.125:6 ... ]
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
title!("Convergence")
plot!(legend=:bottomright)
plot!([extrema(timesteps)... ], repeat([1e-16],2), linestyle=:dash, c="black", alpha=.5, label="")
savefig("$imgdir/rates_error.png")
display(Plots.current())

plot_rates(solvers, timesteps, costs; labels="SOLVER, γ=RATE", ignore=3)
ylabel!("Evaluation time")
title!("Costs")
plot!(legend=:top)
savefig("$imgdir/rates_costs.png")
display(Plots.current())
