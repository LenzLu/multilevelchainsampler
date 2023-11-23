# include("abstract.jl")
# include("dynamics.jl")

using Statistics

## Time evolution

function nodal_deviation(
    grid::PowerGrid, node::Int64,
    perturbation = sample_perturbation(10);
    u_fix = synchronous_state(grid),
    t_final = 180.0, solver=RK4(), kwargs... )

    N = nv(grid.G)
    n_sample = size(perturbation, 1)
    tspan = (0.0,t_final)

    # thread parallel ODE solves
    u_final = zeros(n_sample, length(u_fix))
    #Threads.@threads
    for k=1:n_sample
        u = u_fix[:]
        u[node]   += perturbation[k, 1]
        u[node+N] += perturbation[k, 2]

        problem = ODEProblem(swing_dynamics!, u, tspan, grid)
        solution = solve(problem, solver; kwargs...)
        u_final[k,:] .= solution.u[end];
    end

    ϕ,ω = u_final[:,node], u_final[:,node+N];
    deviation = abs.(ϕ .- u_fix[node])
    return deviation
end

function nodal_basin_stability(
    grid::PowerGrid, node::Int64;
    u_fix = synchronous_state(grid),
    N=100, threshold = 0.05, t_final = 20.0)

    perturbation = sample_qmc_perturbation(N)
    deviation = nodal_deviation(grid, node, perturbation;
                    u_fix, t_final, adaptive=true, abstol=0.5, reltol=0.5)

    return mean( deviation .< threshold )
end

function basin_stability(grid::PowerGrid;
    N=100, threshold = 0.05, t_final = 20.0)

    u_fix = synchronous_state(grid)

    S = zeros(nv(grid.G))
    for i=1:nv(grid.G)
        S[i] = nodal_basin_stability(grid, i)
    end
    return mean( S )
end
