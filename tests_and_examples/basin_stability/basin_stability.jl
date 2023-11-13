include("../../ChainSamplers/all.jl")
include("../../NetworkEnsembles/all.jl")
include("../../PowergridStability/all.jl")

N = 10; p = 0.1  # background ensemble
τ = 2            # proposal change size
K = 6.0; α = 1.0 # dynamics paramters

netsampler = ErdoesRenyiSampler(N, p)
gridsampler = EnergyGridEnsemble(netsampler,τ,K,α)

chain_length = 100
sampler = MetropolisHastings(x -> stability(x; n_sample=10), gridsampler, chain_length)

println("Sampling chain ")
@time X,E = sample_chain(sampler)
println("\n")

using Plots
using GraphPlot

imgdir = "$(@__DIR__)/imgs"
histogram(E, normalized=true)
title!("Energy histogram")
savefig("$imgdir/mh_Ehist.png")
display(Plots.current())

imin = argmin(E)
g = X[imin]
plot_steady(g)
title!("Minimal energy $(E[imin])")
savefig("$imgdir/mh_Emin.png")
display(Plots.current())


imax = argmax(E)
g = X[imax]
plot_steady(g)
title!("Maximal energy $(E[imax])")
savefig("$imgdir/mh_Emax.png")
display(Plots.current())


#=
L = 2
samplesizes = 10 .+ 2 .^ [ L:-1:1 ... ]
surrogates = [
    x -> stability(x; n_sample=n)
    for n in samplesizes
]
sublengths = repeat([10], L-1)

chain_length = 100
energy = surrogates[1]; surrogates = Function[ surrogates[2:end] ... ]
sampler = DelayedAcceptanceMetropolisHastings(energy, surrogates, gridsampler, chain_length, sublengths)
X,E = sample_chain(sampler)
=#
