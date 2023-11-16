using Graphs
using Statistics
using Plots
using GraphRecipes

using Revise
using MultilevelChainSampler

function aspl(g::Graph)
    mean([ mean(dijkstra_shortest_paths(g, i).dists) for i=1:nv(g) ])
end
run(`clear`); println(repeat("=",80)... )
println("Sampling by average shortest path length (aspl.)")

imgdir = "$(@__DIR__)/imgs"
for N = [5,10,50]
    R = ErdoesRenyiSampler(N, 2/(N-1))
    S = MetropolisHastings(aspl, R, 100000)

    println("N = $N nodes")
    x0 = path_graph(N)
    t = @elapsed X,E = sample_chain(S; x0)
    println("Sampling took $t seconds.")

    histogram(E, normalized=true)
    title!("N=$N, p=$(2/(N-1))aspl. histogram")
    display(Plots.current())

    indices = sortperm(E)
    stride = S.length ÷ 20; println(stride)
    indices = [ indices[1:stride:end-1]... , indices[end] ]
    for (i, idx) in enumerate(indices)
        x,e = X[idx],E[idx]
        graphplot(x, curves=false)
        title!("(N=$N) aspl. <l>=$(round(e,digits=5)) \n")
        savefig("$imgdir/aspl_N$(string(N,pad=2))_$(string(i,pad=3)).png")
        display(Plots.current())
    end
end
