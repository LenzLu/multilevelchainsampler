using LinearAlgebra
using OrdinaryDiffEq
using Graphs
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



g = PowerGrid(path_graph(2), [-1,1], 10.0, 0.1)
u_fix = synchronous_state(g); threshold=0.05; t_final = 20.0

N = length(u_fix)÷2
v₀,F = MultilevelChainSampler.sensitivity_inplace(u_fix, swing_dynamics!, MultilevelChainSampler.swing_dynamics_derivative!)

problem = ODEProblem(F, v₀, (0, t_final), g)
println("")
timed = @elapsed solution_fine  = solve(problem, RK4(); adaptive=true, abstol=1e-9, reltol=1e-9)
println("Fine solution ", timed)
timed = @elapsed solution = solve(problem, RK4(); adaptive=true, abstol=0.1, reltol=0.1)
println("Course solution ", timed)

G_fine = solution_fine[end][2N+1:end]
G = solution.u[end][2N+1:end]
error = sum( (G .- G_fine).^2 )

S = basin_stability(g; N=500)
S_proxy = stability_sensitivity_proxy(g)
@test abs(S - S_proxy) < 0.1
