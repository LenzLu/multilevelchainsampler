#=
using Statistics
using IterTools
using OrdinaryDiffEq
using QuasiMonteCarlo
=#

function deviation(g::PowerGrid, perturbation::Vector; t_final=120.0)
    u0 = g.syncstate .+ perturbation
    N = length(u0)÷2 
    problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
    solution = solve(problem, RK4(); adaptive=true, abstol=0.1, reltol=0.1)
    u_final = solution.u[end]

    return ( u_final[N+1:end] .- g.syncstate[N+1:end] )
end

function stable(g::PowerGrid, perturbation::Vector; t_final=120.0, tolerance=0.1)
    δu = deviation(g, perturbation; t_final)
    S = maximum(abs.(δu)) < tolerance
    return S
end

function basin_stability(g::PowerGrid, J::Vector{Int64}, δx=sample_perturbations(length(J), 100); t_final=120.0, tolerance=0.1)
    N = nv(g.grid)
    perturbations = perturb_nodes(N, J, δx)
    perturbations = [perturbations[i, :] for i=axes(perturbations,1)]
    s = δx -> stable(g, δx; t_final, tolerance)
    S = s.(perturbations)
    return S
end

function basin_stability(g::PowerGrid, m::Int64, n::Int64=10; t_final=120.0, tolerance=0.1)
    N = nv(g.grid); b = binomial(N, m)
    S = zeros(b)

    δx=sample_perturbations(length(J), n)
    for (k,J) = enumerate(subsets(1:N, m))
        S[k] .= basin_stability(g, J, δx; t_final, tolerance)
    end
    return mean(S)
end
