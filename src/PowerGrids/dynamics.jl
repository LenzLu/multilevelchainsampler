#=
using LinearAlgebra
using Statistics
using SimpleGraphs
using NonlinearSolve
using OrdinaryDiffEq
=#

function kuramoto_coupling(ϕ::Vector, G::Graph, K::Vector)
    N = length(ϕ)
    F = zeros(eltype(ϕ), N)
    for (i,e) in enumerate(edges(G))
        u,v = e.src, e.dst
        F[u] += K[i] * sin(ϕ[v] - ϕ[u])
        F[v] += K[i] * sin(ϕ[u] - ϕ[v])
    end
    return F
end

function swing_dynamics!(dₜu, u, p, t)
    N = length(u) ÷ 2
    ϕ = u[1:N]; ω = u[N+1:2N]

    G,W,K,α = p.grid, p.power, p.coupling, p.damping 

    dₜu[1:N]    .= ω
    dₜu[1+N:2N] .= W .- α .* ω .+ kuramoto_coupling(ϕ, G, K);
end


function kuramoto_coupling_derivative(ϕ::Vector, G::Graph, K::Vector)
    N = length(ϕ)
    J = zeros(eltype(ϕ), N, N)
    for (i,e) in enumerate(edges(G))
        u,v = e.src, e.dst
        J[u,u] += - K[i] * cos(ϕ[v] - ϕ[u])
        J[u,v] +=   K[i] * cos(ϕ[v] - ϕ[u])
        J[v,u] +=   K[i] * cos(ϕ[u] - ϕ[v])
        J[v,v] += - K[i] * cos(ϕ[u] - ϕ[v])
    end
    return J
end

function swing_dynamics_derivative!(J, u, p, t)
    N = length(u) ÷ 2
    ϕ = u[1:N]; ω = u[N+1:2N]

    G,W,K,α = p.grid, p.power, p.coupling, p.damping 

    J[   1:N,    1:N] .= zeros(N,N)
    J[   1:N, N+1:2N] .= I(N)
    J[N+1:2N,    1:N] .= kuramoto_coupling_derivative(ϕ, G, K)
    J[N+1:2N, N+1:2N] .= - α .* I(N);
end

function is_stable(syncstate::Vector, grid::Graph, power::Vector, coupling::Vector, damping::Vector)
    N = nv(grid); p = (; grid, power, coupling, damping)
    J = zeros(2N,2N)
    swing_dynamics_derivative!(J,syncstate,p,0.0)
    λ = real.( eigen(J[N+1:2N,N+1:2N]).values )
    return maximum(λ) <= 0.0
end

# Steady state
function synchronous_state(grid::Graph, power::Vector, coupling::Vector, damping::Vector)
    N = nv(grid)
    p = (; grid, power, coupling, damping)

    # Nonlinear dynamics solver
    function f(u,p)
        dₜu = zeros(eltype(u), length(u))
        swing_dynamics!(dₜu, u, p, 0)
        dₜu[1:N] .-= mean(dₜu[1:N])
        return dₜu[2:end]
    end

    u = zeros(2N)
    for i=1:10
        problem = NonlinearProblem(f, u, p)
        solution = solve(problem, NewtonRaphson())
        u = solution.u
        if is_stable(u, p...) return u 
        else u = randn(length(u)) end
    end
    @warn "Synchronous state not found!"

end
