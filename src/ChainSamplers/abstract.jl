
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
