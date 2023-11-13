include("abstract.jl")

using LinearAlgebra
using Graphs
using NonlinearSolve
using OrdinaryDiffEq

## Dynamics

function KuramotoCoupling(ϕ::Vector, G::Graph)
    N = nv(G)
    coupling = zeros(eltype(ϕ), N)
    for i=1:N
        for j=neighbors(G, i)
            coupling[i] += sin(ϕ[i] - ϕ[j])
        end
    end
    return coupling
end


function swing_dynamics!(dₜu, u, p::PowerGrid ,t)
    N = length(u)÷2 # Unpack state
    ϕ = u[1:N]; ω = u[N+1:2N]

    # Swing equation
    dₜϕ = ω
    dₜω = p.W .- p.α .* ω .- p.K .* KuramotoCoupling(ϕ, p.G)

    dₜu[1:N] = dₜϕ; dₜu[N+1:2N] =  dₜω  # pack vector
    return dₜu
end



## Solver

function solve_system(grid::PowerGrid,  u₀; t_final=10.0, tol=1e-6)
    problem = ODEProblem(swing_dynamics!, u₀, (0.0, t_final), grid)
    solution = solve(problem, Rodas5P(); abstol=tol, reltol=tol)

    #show_dynamics(solution)
    u = solution.u[end]
end


#= Analytical Derivative of KuramotoCoupling
# function KuramotoCouplingDerivative(ϕ::Vector, G::Graph)
#     N = nv(G)
#     coupling = zeros(eltype(ϕ), N, N)
#     for i=1:N
#         for j=neighbors(G, i)
#             coupling[i,i] += cos(ϕ[i] - ϕ[j])
#             coupling[i,j] -= cos(ϕ[i] - ϕ[j])
#         end
#     end
#     return coupling
# end
# # only used for non-linear angle solver =#


function steady_state(grid::PowerGrid)

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

    # Nonlinear dynamics solver
    function f(u,p)
        dₜu = zeros(eltype(u), length(u))
        swing_dynamics!(dₜu, u, p, 0)
        dₜu[1:N] .-= mean(dₜu[1:N])
        return dₜu[2:end]
    end
    N = nv(grid.G); u = zeros(2N)
    problem = NonlinearProblem(f, u, grid)
    solution = solve(problem, NewtonRaphson())
    u = solution.u

    return u
end
