
function stable(g::PowerGrid, perturbation::Vector; t_final=120.0, tolerance=0.1)
    u0 = g.syncstate .+ perturbation
    N = length(u0)÷2 
    problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
    
    max_deviation = 3 * max( tolerance, maximum(abs.(perturbation)) ) 
    desynced(u,t,integrator) = ( maximum(abs.(u[N+1:2N] .- g.syncstate[N+1:2N])) > max_deviation )
    callback = ContinuousCallback(desynced, terminate!)
    
    solution = solve(problem, RK4(); adaptive=true, abstol=1e-2, reltol=1e-1, callback)
    u_final = solution.u[end]

    δu = ( u_final[N+1:end] .- g.syncstate[N+1:end] )
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
        S_J = basin_stability(g, J, δx; t_final, tolerance)
        S[k] = mean(S_J)
    end
    return mean(S)
end
