include("../../ChainSamplers/all.jl")
include("../../NetworkEnsembles/all.jl")


using Distributions
function energy(g::Graph)
    c = betweenness_centrality(g)
    e = 1 .- mean(c)
    return e
end

R = ErdoesRenyiSampler(20, 0.1)
S = MetropolisHastings(energy, R, 10000)
X,E = sample_chain(S)

using Plots
histogram(E)
display(Plots.current())

using GraphPlot

imin = argmin(E)
plt = gplot(X[imin])
#title!("Minimal energy graph")
display(plt)

imax = argmax(E)
plt = gplot(X[imax])
#title!("Maximal energy graph")
display(plt)
