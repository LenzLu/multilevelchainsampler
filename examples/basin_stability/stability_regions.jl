using Graphs
using OrdinaryDiffEq
using Statistics
using Plots
using GraphRecipes

using Revise
using MultilevelChainSampler
imgdir = "$(@__DIR__)/imgs"

function plot_grid_stability(g::PowerGrid,
    c_passed = RGBA(1, .5, 0, .8); c_failed = RGBA(0, 0, .6, .2))

    perturbation = grid_perturbation(18, 18)
    #perturbation = sample_qmc_perturbation(324) # QMC sampling
    deviation = MultilevelChainSampler.nodal_deviation(g,1,perturbation;adaptive=true, abstol=0.01,reltol=0.01)
    idx = deviation .< 0.05; nidx = .! idx

    plt_stability = plot()
    #scatter!(perturbation[:,1], perturbation[:,2], color="gray", alpha=.2, label="")
    scatter!(perturbation[idx,1], perturbation[idx,2], color=c_passed, label="passed")
    scatter!(perturbation[nidx,1], perturbation[nidx,2], color=c_failed, label="failed")
    xlabel!("Δϕ₁"); ylabel!("Δω₁")
    S = round( mean(idx), sigdigits=3 )
    title!("Stability S=$S")

    annotate!(6, -2., "K = $(round(g.K, sigdigits=2)), α=$(round(g.α, sigdigits=2))", fontsize=1)
    plt_topology = plot(g.G, names=1:nv(g.G), fontsize=10, nodesize=.15, curves=false)
    plot(plt_stability, plt_topology)
    return Plots.current()

end


K = 2.0; a = 0.5
g = PowerGrid(complete_graph(2), [1,-1], K, a)

plot_grid_stability(g)
savefig("$imgdir/stability_K2.png")
display(Plots.current())



N=19; K = 9.0; a = 0.5
ensemble = PowerGridEnsemble(ErdoesRenyiSampler(N,0.15),10,K,a)
g = initialize(ensemble)
plot_grid_stability(g)

propose!(ensemble, g)
plot_grid_stability(g)

timed_single = @elapsed nodal_basin_stability(g, 1)
timed_total = @elapsed basin_stability(g)

println("Single run ", timed_single, " × ", nv(g.G), " = ", timed_single*nv(g.G))
println("Total run ", timed_total)
