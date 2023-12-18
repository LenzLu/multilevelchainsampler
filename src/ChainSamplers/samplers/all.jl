include("chain.jl")

# Chain sampling algorithm
abstract type AbstractChainSampler end
function sample(e::AbstractEnergyFunction, s::AbstractChainSampler, g::AbstractProposalGenerator; x0=propose(g))
    E0 = evaluate(e, x0)
    return ChainSample([x0], [E0], [0], [nothing])
end
function sample(f::Function, s::AbstractChainSampler, g::AbstractProposalGenerator; kwargs...)
    e = SimpleEnergyFunction(f)
    return sample(e, s, g; kwargs...)
end

include("flathistogram.jl")
include("metropolis.jl")
