#using Plots
include("dynamics.jl")

using Revise
using MultilevelChainSampler

# Naive estimator
dt = 1e-2
δu = sample_perturbation(100)
q = quantity_of_interest(δu, dt)
μ,σ = mean(q), std(q)
print("Naive $μ ± $σ")


# Fixed size estimator
μ,σ = multilevel_estimator(
    sample_perturbation, quantity_of_interest,
    [512, 256, 128, 64, 32, 16, 8, 4, 2]; h0=1.0, m=0.5)
print("Fixed sizes $μ ± $σ")

#=
μ,σ = adaptive_multilevel_estimator(
    sample_perturbation, quantity_of_interest;
    h0=1.0, m=0.5, α=?, γ=?)
=#
