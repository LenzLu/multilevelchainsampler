using Makie, CairoMakie
using Statistics
using Random 
using Revise
using MultilevelChainSampler

Random.seed!(42)
GaussianEnergy(n::Vector{Int64}) = TrackedEnergyFunction(
    SamplingBasedEneryFunction(
        rand(0.1..9.9, 1000), n, (x,s) -> (x - s)^2
    )
)

function plot_energy!(ax,f; kwargs...)
    t = collect(0:.01:10)
    f = exp.(-f.(t))
    f = f ./ ( mean(f) * 10) 
    lines!(ax, t, f; kwargs...)
end

function plot_chain_statistics!(fig, chain; unfolded=true, args...)
    #=
    ax = Axis(fig[1,3])
    hist!(ax, chain[:acceptance], normalization=:pdf, args...)
    =#
    
    if unfolded   
        chain = unfold(chain)
    end
    ax = Axis(fig[1,1]); ax.xlabel = L"State $x$"

    hist!(ax, chain.states, normalization=:pdf, args...)
    ax_state = ax

    ax = Axis(fig[1,2]); ax.xlabel = L"Energy $E$"
    
    hist!(ax, chain.energies, normalization=:pdf, args...)
    xlims!(extrema(chain.energies)...)
    ax_energy = ax

    return ax_state, ax_energy
end
