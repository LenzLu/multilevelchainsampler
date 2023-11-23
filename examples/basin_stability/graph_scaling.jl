using Graphs
using OrdinaryDiffEq
using Statistics
using Plots


using Revise
using MultilevelChainSampler

clearconsole(); println(repeat("=",80)... )
println("Measuring graph scaling")
imgdir = "$(@__DIR__)/imgs"


K = 5.0; a = 0.5
TIMER_ITERS = 10

nodes = [10:2:16...,  18:8:64 ..., 66:16:128 ... ]
tolerances = [1e-2, 1e-3, 1e-4]
rk_costs  = zeros(length(nodes), length(tolerances))
rk_errors = zeros(length(nodes), length(tolerances))
#=stepsizes = [0.075, 0.01]
ie_costs  = zeros(length(nodes), length(stepsizes))
ie_errors  = zeros(length(nodes), length(stepsizes)) =#
timespan = (0.0, 10.0)


for (i,N) in enumerate(nodes)

    println("N = $N")

    G = complete_graph(N)
    for v = N:-1:convert( Int64, ceil(N*.5))+1
        rem_edge!(G, 1, v)
    end

    W = ones(nv(G)); W[1:nv(G)÷2] .= -1
    g = PowerGrid(G,W,K,a)
    u_fix = synchronous_state(g)
    perturbation = [ π/2, 10.0 ]
    u₀ = u_fix[:]; u₀[ [1,N+1] ] += perturbation
    problem = ODEProblem(swing_dynamics!, u₀, timespan, g)

    print("Evaluating fine solution")
    @time fine_solution = solve(problem, RK4(); adaptive=true, reltol=1e-12)
    u_fine = fine_solution.u[end][1]

    for (j,tol) in enumerate(tolerances)
        solver = RK4(); kwargs = (; adaptive=true, abstol=tol, reltol=tol )
        solution = solve(problem, solver; kwargs...)
        u = solution.u[end][1]

        rk_errors[i,j] = abs(u_fine - u) / abs(u_fine)
        rk_costs[i,j] = mean([ @elapsed solve(problem, solver; kwargs...) for k=1:TIMER_ITERS ])
    end
    #=
    for (j,h) in enumerate(stepsizes)
        solver = ImplicitEuler(); kwargs = (; adaptive=false, dt=h )
        #solver = ImplicitEuler(); kwargs = (; adaptive=true, abstol=h, reltol=h )
        solution = solve(problem, solver; kwargs...)
        u = solution.u[end][1]
        ie_errors[i,j] = abs(u_fine - u) / abs(u_fine)
        ie_costs[i,j] = mean([ @elapsed solve(problem, solver; kwargs...) for k=1:TIMER_ITERS ])
     end
     =#
end

plot()
for (j,tol) in enumerate(tolerances)
    x = j / length(tolerances); c = RGB(.5,.5 * x, .5 * (1-x))
    plot!(nodes, rk_costs[:,j], label="RK4, TOL=$(round(tol, sigdigits=3))"; c)
end
#=
for (j,h) in enumerate(stepsizes)
    x = j / length(stepsizes); c = RGB(0,.5 + .5 * x, .5 + .5 * (1-x))
    plot!(nodes, ie_costs[:,j], label="ImplicitEuler, h=$(round(h, sigdigits=3))"; c)
end
=#

plot!(nodes, 1e-5 .* nodes .^ 2, color="black", linestyle=:dash, label="\$O(N^2)\$")
xlabel!("Number of nodes \$N\$"); ylabel!("Computation time")
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
plot!(legend=:bottomright)
title!("Costs per network size")
savefig("$imgdir/scaling_nv_vs_cost.png")
display(Plots.current())

plot()
for (j,tol) in enumerate(tolerances)
    x = j / length(tolerances); c = RGB(.5,.5 * x, .5 * (1-x))
    plot!(nodes, rk_errors[:,j], label="RK4, TOL=$(round(tol, sigdigits=3))"; c)
end
#=
for (j,h) in enumerate(stepsizes)
    x = j / length(stepsizes); c = RGB(0,.5 + .5 * x, .5 + .5 * (1-x))
    plot!(nodes, ie_errors[:,j], label="ImplicitEuler, h=$(round(h, sigdigits=3))"; c)
end
=#
plot!(nodes, 1e-5 .* nodes .^ 2, color="black", linestyle=:dash, label="\$O(N^2)\$")
#plot!(nodes, 1e-5 .* nodes .^ 1, color="black", linestyle=:dash, label="\$O(N)\$")
xlabel!("Number of nodes \$N\$"); ylabel!("Relative error")
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
plot!(legend=:bottomright)
title!("Error per network size")
savefig("$imgdir/scaling_nv_vs_error.png")
display(Plots.current())
