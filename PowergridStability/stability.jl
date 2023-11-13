include("abstract.jl")
include("dynamics.jl")

using QuasiMonteCarlo
using Statistics


function stability(grid::EnergyGrid, u_steady::Vector, node::Integer; n_sample=100, threshold=0.05)

    ## Sample random perturbations
    #U = rand(n_sample)
    U = QuasiMonteCarlo.sample(n_sample, 1, SobolSample())
    Δϕ = ( 2 .* U .- 1) .* π * 0.1

    Δω = randn(n_sample) .* 0.0

    N = nv(grid.G)
    ϕ_fix = u_steady[1:N]
    failures = 0
    for i=1:n_sample
        u = copy(u_steady)
        u[node]   = mod(u[node] + Δϕ[i], 2π)
        u[node+N] = u[node+N] + Δω[i]
        u = solve_system(grid, u)
        ϕ = u[1:N]; ω = u[N+1:2N]

        if maximum(abs.(ϕ - ϕ_fix)) >= threshold
            failures += 1
        end
    end
    return 1.0 - failures / n_sample
end

function stability(grid::EnergyGrid; n_sample=10)
    N = nv(grid.G)
    u = steady_state(grid)

    S = [stability(grid, u, i; n_sample) for i=1:N]
    #return
    mean(S)
end

#===== Test ====================================================================
run(`clear`)
N = 20; K = 6.0; α = 0.1
net = ErdoesRenyiSampler(10, 0.1)
S = EnergyGridEnsemble(net, K, α)
grid = rand(S)

#ϕ = zeros(N)
#ω = zeros(N)
#u = steady_state(grid)
s = stability(grid)
print(s)
===============================================================================#
