using CairoMakie
using Graphs
using Random

using Revise
using MultilevelChainSampler

datadir = "$(@__DIR__)/data"

function plot_chain!(ax, s::ChainSample; unfolded=true)
    if unfolded 
        s = unfold(s)
    end
    hist!(ax, s.energies)
end

@time s = load_chain("$datadir/chain_mh.jld")
fig = Figure(); ax = Axis(fig[1,1])
plot_chain!(ax, s); ax.title = "MetropolisHastings"
display(fig)

@time s = load_chain("$datadir/chain_wl.jld")
fig = Figure(); ax = Axis(fig[1,1])
plot_chain!(ax, s); ax.title = "WangLandau"
display(fig)
