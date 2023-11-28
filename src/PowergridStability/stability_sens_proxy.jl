using OrdinaryDiffEq



function stability_sensitivity_proxy(g::PowerGrid;
    u_fix = synchronous_state(g), threshold=0.05, t_final = 20.0)

    N = length(u_fix)÷2
    v₀,F = sensitivity_inplace(u_fix, swing_dynamics!, swing_dynamics_derivative!)

    problem = ODEProblem(F, v₀, (0, t_final), g)
    solution = solve(problem, RK4(); adaptive=true, reltol=1e-3, abstol=1e-3)

    #@assert sum( (solution.u[end][1:2N] .- u_fix ) .^2 ) < 0.1 "No synchronous state found!"
    G = reshape( solution.u[end][2N+1:end], 2N, 2N)

    S = zeros(N)
    for i=1:N
        S[i] = min(π, threshold / sqrt( sum( G[:,i] .^ 2 ) )) / π
    end
    return mean(S)
end
