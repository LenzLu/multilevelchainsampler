using Distributions
using QuasiMonteCarlo

## Stochastic Perturbation

function sample_perturbation(N=10)
    phase = Uniform(-π,π)
    velocity = Normal(0.0, 2.0)
    δϕ = rand(phase, N)
    δω = rand(velocity, N)
    perturbation = hcat( δϕ, δω )
end

function sample_qmc_perturbation(N=10)
    U = QuasiMonteCarlo.sample(N,2, SobolSample())
    δϕ = π .* ( 2 .* U[1,:] .- 1)
    δω = 2.0 .* ( 2 .* U[2,:] .- 1)
    perturbation = hcat( δϕ, δω )
end
## Determinisitic Perturbation

function grid_perturbation(angles=48, velocities=5)
    U = [1:angles ... ] ./ angles
    V = [1:velocities ... ] ./ velocities

    δϕ = @. (2.0 * U - 1.0) * π;
    δω = @. (2.0 * V - 1.0) * 2.0;
    δϕ = repeat(δϕ, outer=velocities)
    δω = repeat(δω, inner=angles)

    perturbation = hcat(δϕ, δω)
end
