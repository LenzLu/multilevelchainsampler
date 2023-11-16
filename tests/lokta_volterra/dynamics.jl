using OrdinaryDiffEq
using Statistics
#using NonlinearSolve

using Revise
using MultilevelChainSampler

function lokta_volterra!(dₜu, u, p, t)
    a,b,c,d = p
    dₜu[1] =  a * u[1]  - b * u[1]*u[2]
    dₜu[2] = -c * u[2]  + d * u[1]*u[2]
end

p = [1,1,1,1]
u₀ = [.5, 2]
δ  = 0.2
t_final  = 6.0

function sample_perturbation(n::Int64)
    δu = (2 .* rand(n, 2) .- 1) .* δ
    return δu
end

function deviation(δu; solver=Euler(), adaptive=false, dt=0.1, kwargs... )
    problem = ODEProblem(dynamics!, u₀ .+ δu, (0, t_final), p)
    solution = solve( problem, solver; adaptive, dt, kwargs... )
    hcat(solution.u ... )
end


#function deviation()
#=
function quantity(δu, dt; solver = Euler(), adaptive=false)
    nsamples = size(δu,1)
    u = zeros(nsamples)

    for i=1:nsamples
        problem = ODEProblem(lokta_volterra!, u₀ .+ δu[i,:], (0,t_final), p)
        solution = solve(problem; solver, adaptive, dt)
        u[i] = solution[end][1]
    end
    return u
end
=#
