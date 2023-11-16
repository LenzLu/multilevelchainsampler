using Distributions

## Stochastic Perturbation

function sample_perturbation(nsample=10)
    phase = Uniform(-π,π)
    velocity = Normal()

    δϕ = rand(phase, nsample)
    δω = rand(velocity, nsample)
    perturbation = hcat( δϕ, δω )
end

## Determinisitic Perturbation

function grid_perturbation(angles=48, velocities=5)
    U = [0:angles-1 ... ] ./ angles
    V = [0:velocities-1 ... ] ./ velocities

    δϕ = @. (2*U - 1) * π;
    δω = @. (2*V - 1)
    δϕ = repeat(δϕ, inner=1, outer=velocities)
    δω = repeat(δω, outer=angles, inner=1)

    perturbation = hcat(δϕ, δω)
end
