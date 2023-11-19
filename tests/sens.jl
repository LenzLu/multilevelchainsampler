using OrdinaryDiffEq
using Test

using Revise
using MultilevelChainSampler


f(u,p,t) = [-u[2],u[1]]; u₀ = [1.0, 0.0]
∂ᵤf(u,p,t) = [ 0 -1 ; 1 0 ]
# ⟹ u(t) = [cos(t), sin(t)]


v₀,F = sensitivity(u₀, f, ∂ᵤf)
@test size(v₀) == (6,)
@test size(F(v₀,[],0)) == (6,)

problem = ODEProblem(f, u₀, (0.0, 3.1415))
solution = solve(problem, Tsit5())
u_final = solution.u[end]
@test sum(abs.( u_final .+ u₀ )) < 1e-3

problem = ODEProblem(F, v₀, (0, 3.1415))
solution = solve(problem, Tsit5())
v_final = solution.u[end]
@test maximum(abs.(v_final[3:end] .+ [I(2)... ])) < 1e-3
