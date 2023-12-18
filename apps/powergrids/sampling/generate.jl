using Statistics
using Graphs
using Random

using Revise
using MultilevelChainSampler

datadir = "$(@__DIR__)/data"

q = ErdoesRenyiEnsemble(20,0.12)
q = PowerGridProposal(q)
e = EnergyFunction(r_uni)

sampler_mh = MetropolisHastings(100000)
@time s_mh = sample(e, sampler_mh, q)
@time save_chain(s_mh, "$datadir/chain_runi_mh.jld")

sampler_wl = WangLandau(0.95, 1.0, 100000; nbins=200)
@time s_wl = sample(e, sampler_wl, q)
@time save_chain(s_wl, "$datadir/chain_runi_wl.jld")
