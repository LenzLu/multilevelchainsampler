include("../dynamics.jl")

using Plots
using StatsPlots
using Revise
using MultilevelChainSampler
imgdir = "$(@__DIR__)/../imgs"


blend(c1, c2, t) = c1 .* (1-t) .+ c2 .* t

# Distribution Histograms
x = sample_perturbation(5000)
stephist(u₀[1] .+ x[:,1], alpha=.25, normalize=true, width=1, color=:black, label="Initial")
timesteps = [0.5, 0.33, 0.25, 0.1, 0.05, 0.01]
for (i,dt) = enumerate(timesteps)
    q = quantity_of_interest(x, dt)
    n = length(timesteps); color = RGB(i/n, 1-i/n, 1-i/n)
    stephist!(q, alpha=.2, normalized=true, label="dt=$dt", width=4; color)
end
title!("Distribution mapping")
savefig("$imgdir/distributions.png")
display(Plots.current())

# Statistics
plot()
timesteps = 0.1 .^ [0:.1:3 ... ]
xlims!(extrema(timesteps)... )
nsamples = [2,10]; K = 10
for (i,N) = enumerate(nsamples)
    m = zeros(length(timesteps), K)
    s = zeros(length(timesteps), K)
    for k=1:K
        for (j,dt) in enumerate(timesteps)
            x = sample_perturbation(N)
            q = quantity_of_interest(x, dt)
            m[j,k] = mean(q)
            s[j,k] = std(q)
        end
    end
    n = length(nsamples)
    cm = RGB(blend([1,0,0], [.5,0,0], i/n)...)
    cs = RGB(blend([0,0,1], [0,0,.5], i/n)...)
    errorline!(timesteps, m, color=cm, label="mean N=$N")
    errorline!(timesteps, s[:,:], color=cm, label="std N=$N")
end
xlabel!("Time step Δt")
title!("Statistics")
#xaxis!(:log); yaxis!(:log)
savefig("$imgdir/statistics.png")
display(Plots.current())
