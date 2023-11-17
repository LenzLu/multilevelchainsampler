## Proposal generator

abstract type ProposalGenerator{StateType} end

function initialize(::ProposalGenerator{T})::T where{T} error("Not implemented") end
function propose!(::ProposalGenerator{T}, ::T) where{T} error("Not implemented") end

# Default transition(s)
function transition(::ProposalGenerator{T}, ::T, ::T) where{T} 1.0 end
function transitions(g::ProposalGenerator{T}, x::T, y::T) where{T} return (transition(g,x,y), transition(g,y,x)) end


## Chain Sampler

abstract type ChainSampler{StateType} end

function sample_chain(::ChainSampler{T}) where {T} error("Not implemented") end



const FunctionIterable = Union{Vector{<:Function}, Tuple{Vararg{Function}}}


include("metropolis_hastings.jl")
include("mh_delayed_acceptance.jl")
include("mh_multilevel.jl")
#include("wanglandau.jl")
#include("histutils.jl")
