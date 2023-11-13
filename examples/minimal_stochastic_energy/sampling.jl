include("../../src/all.jl")
using Statistics

## Trivial Generator

mutable struct Dummy
    x::Float64 # State
end
struct RandomWalk <: ProposalGenerator{Dummy}
    s::Float64 # Step size
end
initialize(R::RandomWalk) = Dummy(randn() * R.s)
propose!(R::RandomWalk, d::Dummy) = ( d.x += randn() * R.s; d)

## Gaussian target distribution

function energy(d::Dummy) return d.x^2 end
function energy_surrogate(d::Dummy; c=1.0, n=1)
    z = c.* ( rand(n) .- 0.5 )
    return mean( (d.x .- z) .^ 2 )
end


using Plots
function create_plots(sampler::ChainSampler, name::String)
    imgdir = "$(@__DIR__)/imgs"

    println("Sampling $name :");

    @time chain_X, chain_E = sample_chain(sampler, x0=Dummy(10))
    print("\n")


    N = length(chain_X)
    t = [1:N ...]
    chain_X = [d.x for d in chain_X]

    # Trajectories

    plot(t, chain_X)
    title!("$name State Trajectory")
    savefig(Plots.current(), "$(imgdir)/StateTrajectory_$(name).png")
    display(Plots.current())
    plot(t, chain_E)
    title!("$name Energy Trajectory")
    savefig(Plots.current(), "$(imgdir)/EnergyTrajectory_$(name).png")
    display(Plots.current())


    # Sample statistics
    m = [ mean(chain_X[1:n])             for n=1:N ]
    s = [ mean((chain_X[1:n] .- m[n]).^2) for n=1:N ]
    s = sqrt.(s)

    plot(t, m, label="running mean", c=:red)
    plot!([t[1], t[end]], repeat([m[end]],2), label="final mean $(m[end])", c=:red , ls=:dash)
    plot!(t, s, label="running std", c=:cyan)
    plot!([t[1], t[end]], repeat([s[end]],2), label="final std $(s[end])", c=:cyan, ls=:dash)
    title!("$name Statistics")
    savefig(Plots.current(), "$(imgdir)/Statistics_$(name).png")
    display(Plots.current())


    # Compare to target distribution
    histogram(chain_X, normalize=true, label="state histogram")
    # Target probability distribution
    h = 0.01
    x = [ -8:h:8 ... ]
    e = energy.(Dummy.(x))
    p = exp.(-e)
    ∫p = h * sum(p)
    p = p / ∫p
    plot!(x, p, label="target distribution", lw=2.5)
    title!(
        "$name histogram"
        #* "where n \$p(x) = \\frac{\\exp(-E(x))}{\\int \\exp(-E(\\xi)) d\\xi}\$"
    )
    savefig(Plots.current(), "$(imgdir)/Target_$(name).png")
    display(Plots.current())


end

# = TEST

run(`clear`); println(repeat("=",80)... )
println("Minimal Stochastical tests ")

R = RandomWalk(0.05)

chain_length = 10000
sampler = MetropolisHastings(energy, R, chain_length)
create_plots(sampler, "Vanilla Metropolis Hastings")

chain_length = 10000
surrogates = Function[ x->energy_surrogate(x; n=16) ]
sampler_TLDA = DelayedAcceptanceMetropolisHastings(energy, Function[energy_surrogate], R, chain_length, [4])
create_plots(sampler, "Two-level Delayed Acceptance")

chain_length = 10000
surrogates = Function[ x->energy_surrogate(x; n=k) for k = [8,4] ]
surrlengths = [16,8]
sampler = DelayedAcceptanceMetropolisHastings(energy, surrogates, R, chain_length, surrlengths)
create_plots(sampler, "Three-level Delayed Acceptance")

chain_length = 10000
surrogates = Function[ x->energy_surrogate(x; n=k) for k = [8,4,2] ]
surrlengths = [16,8,4]
sampler = DelayedAcceptanceMetropolisHastings(energy, surrogates, R, chain_length, surrlengths)
create_plots(sampler, "Four-level Delayed Acceptance")



#chain_length = 10000
#surrogates = [ x->energy_surrogate(x; n=k) for k = [2,4] ]
#sampler_MLDA = DelayedAcceptanceMetropolisHastings(energy, surrogates, R, chain_length, surrlengths)
#create_plots(sampler_MLDA, "Four-level Delayed Acceptance")
