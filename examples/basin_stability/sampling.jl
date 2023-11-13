include("../../src/all.jl")

using Plots
using GraphPlot

N = 10; p = 0.1  # background ensemble
τ = 2            # proposal change size
K = 6.0; α = 2.0 # dynamics paramters

netsampler = ErdoesRenyiSampler(N, p)
gridsampler = PowerGridEnsemble(netsampler,τ,K,α)

imgdir = "$(@__DIR__)/imgs"
function create_plots(sampler, fname)
    println("Sampling chain ")
    @time X,E = sample_chain(sampler)
    println("\n")

    histogram(E, normalized=true)
    title!("Energy histogram")
    savefig("$imgdir/$(fname)_Ehist.png")
    display(Plots.current())

    for (i,(x,e)) in enumerate(histsample(X,E))
        plot_steady(g)
        title!("Grid stability $(e)" )
        savefig("$imgdir/$(fname)_Ebin$(i).png")
        display(Plots.current())
    end
end

chain_length = 10
sampler = MetropolisHastings(x -> stability(x), gridsampler, chain_length)
create_plots(sampler, "mh")
