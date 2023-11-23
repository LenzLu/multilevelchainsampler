# include("abstract.jl")

using LinearAlgebra
using Graphs
using NonlinearSolve
using OrdinaryDiffEq

## Dynamics

function kuramoto_coupling(ϕ::Vector, G::Graph)
    N = nv(G)
    coupling = zeros(eltype(ϕ), N)
    for i=1:N
        for j=neighbors(G, i)
            coupling[i] += sin(ϕ[j] - ϕ[i])
        end
    end
    return coupling
end


function swing_dynamics!(dₜu, u, p::PowerGrid ,t)
    N = length(u)÷2 # Unpack state
    ϕ = u[1:N]; ω = u[N+1:2N]

    # Swing equation
    dₜϕ = ω
    dₜω = p.W .- p.α .* ω .+ p.K .* kuramoto_coupling(ϕ, p.G)

    dₜu[1:N] .= dₜϕ; dₜu[N+1:2N] .=  dₜω  # pack vector
    return dₜu
end


#Analytical Derivative of KuramotoCoupling
function kuramoto_coupling_derivative(ϕ::Vector, G::Graph)
     N = nv(G)
     J_coupling = zeros(eltype(ϕ), N, N)
     for i=1:N
         for j=neighbors(G, i)
             J_coupling[i,i] -= cos(ϕ[j] - ϕ[i])
             J_coupling[i,j] += cos(ϕ[j] - ϕ[i])
         end
     end
     return J_coupling
 end

function swing_dynamics_derivative!(J, u, p::PowerGrid ,t)
    N = length(u)÷2 # Unpack state
    ϕ = u[1:N]; ω = u[N+1:2N]

    # Block definition of Jacobian
    J[   1:N,    1:N] .= zeros(N,N)
    J[   1:N, N+1:2N] .= I(N)
    J[N+1:2N,    1:N] .= p.K .* kuramoto_coupling_derivative(ϕ, p.G)
    J[N+1:2N, N+1:2N] .= - p.α .* I(N);

end

# Steady state
function synchronous_state(grid::PowerGrid)

    # # Linearized angle solver
    # L = laplacian_matrix(grid.G)
    # M = pinv(Matrix(L))
    # ϕ = - 1 / grid.K .* M * grid.W
    # ϕ = mod.(ϕ, 2π)
    # ω = zeros(eltype(ϕ), length(ϕ))
    # u = [ϕ..., ω... ]

    # #  Nonlinear angle solver
    # N = nv(grid.G); ϕ = zeros(N)
    # f(ϕ,p) = p.K * KuramotoCoupling(ϕ, p.G) + p.W
    # j(ϕ,p) = p.K * KuramotoCouplingDerivative(ϕ, p.G)
    # func = NonlinearFunction(f, jac=j)
    # problem = NonlinearProblem(func, ϕ, grid)
    # solution = solve(problem, NewtonRaphson())
    # ϕ = solution.u
    # ω = zeros(eltype(ϕ), length(ϕ))
    # u = [ϕ..., ω... ]

    N = nv(grid.G)

    # Nonlinear dynamics solver
    function f(u,p)
        dₜu = zeros(eltype(u), length(u))
        swing_dynamics!(dₜu, u, p, 0)
        dₜu[1:N] .-= mean(dₜu[1:N])
        return dₜu[2:end]
    end

    for i=1:10
        u = zeros(2N)
        problem = NonlinearProblem(f, u, grid)
        solution = solve(problem, NewtonRaphson())
        u = solution.u
        if is_stable(grid, u) return u end
    end
    @warn "Synchronous state not found!"
end

function is_stable(grid::PowerGrid, u_fix::Vector)
    @assert length(u_fix) == 2*nv(grid.G)
    N = nv(grid.G)

    J = zeros(2N,2N)
    swing_dynamics_derivative!(J,u_fix,grid,0.0)
    λ = real.( eigen(J).values )
    return minimum(λ) <= 0.0
end
