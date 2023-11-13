include("../../src/all.jl")

using Distributions
using Plots
using GraphRecipes

function energy(g::Graph)
    aspl = mean([ mean(dijkstra_shortest_paths(g, i).dists) for i=1:nv(g) ])
    e = aspl
    return e
end



run(`clear`); println(repeat("=",80)... )
println("Minimal Network tests ")

R = ErdoesRenyiSampler(20, 0.1)
S = MetropolisHastings(energy, R, 10000)

X,E = sample_chain(S)

histogram(E, normalized=true)
display(Plots.current())

for (x,e) in histsample(X,E; n_bins=15)
    graphplot(x, curves=false)
    title!("energy E=$e \n")
    display(Plots.current())
end
