# include("abstract.jl")
# include("dynamics.jl")

using Statistics
using QuasiMonteCarlo

## Time evolution

function nodal_deviation(
    grid::PowerGrid, node::Int64;
    threshold = 0.1,
    u_fix = synchronous_state(grid),
    perturbation = sample_perturbation(10),
    t_final = 10.0, solver=Tsit5(), kwargs... )

    N = nv(grid.G)
    n_sample = size(perturbation, 1)
    tspan = (0.0,t_final)

    # thread parallel ODE solves
    u_final = zeros(n_sample, length(u_fix))
    Threads.@threads for k=1:n_sample
        u = u_fix[:]
        u[node]   += perturbation[k, 1]
        u[node+N] += perturbation[k, 2]

        problem = ODEProblem(swing_dynamics!, u, tspan, grid)
        solution = solve(problem, solver; kwargs...)
        u_final[k,:] .= solution.u[end];
    end

    ϕ,ω = u_final[node], u_final[node+N];
    deviation = abs(ϕ - u_fix[node])
    #            + abs(ω - u_fix[node+N])
    #deviation .< threshold
    return deviation
end


function nodal_basin_stability(
    grid::PowerGrid, node::Int64;
    threshold = 0.1,
    u_fix = synchronous_state(grid),
    t_final=10.0)

    # stability function dependent on time step
    solver = Tsit5()
    stability(dt) = Δu -> mean(
        nodal_deviation(grid, node;
          threshold, u_fix, perturbation=Δu,
          solver, t_final, adaptive=false, dt)
    )

    # multilevel monte-carlo estimator
    levels = [1e-1, 1e-2, 1e-3, 1e-4]
    nsamples = [100, 50,  25, 10]
    p(nsample) = sample_perturbation(nsample)
    q = [ stability(dt) for dt in levels ]

    S = multilevel_estimator(p, q; nsamples)
    return S
end


## Stability measures


function basin_stability( grid::PowerGrid; threshold=0.1, t_final=10.0)
    u_fix = synchronous_state(grid)
    S = zeros( nv(grid.G) )
    for node=1:nv(grid.G)
        S[node] = basin_stability(grid, node; threshold, u_fix, t_final)
    end
    return mean(S)
end


#=
function warn_runtime(msg, time)
    if evaltime > 1.0
        if evaltime < 60.0
            evaltime = "$evaltime s"
        else
            evaltime = "$(evaltime ÷ 60.0) min"
        end
        @warn "Stability calculation est. $evaltime"
    end
end

function basin_stability(grid::PowerGrid;
    t_final = 10.0,
    u_fix = synchronous_state(grid),
    solver=RK4(), kwargs...)


    # if verbose
    #     u = u_fix .+ randn(length(u_fix))
    #     problem = ODEProblem(swing_dynamics!, u, (0,t_final), grid)
    #     solvetime = @elapsed solve(problem, solver; kwargs... )
    #     evaltime = N * n_sample * solvetime
    #     warn_runtime(evaltime)
    # end


    S = zeros(size(perturbation,1))
    for i=1:nv(grid.G)
        S[i] = nodal_basin_stability(grid, i;
                t_final, u_fix, perturbation, solver, kwargs...)
    end
    return mean(S)
end
=#
