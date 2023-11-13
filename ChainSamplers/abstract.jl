abstract type ProposalGenerator{T} end
function initialize(::ProposalGenerator{T})::T where{T} end
function propose!(::ProposalGenerator{T}, ::T) where{T} end
function probability(::ProposalGenerator{T}, ::T, ::T) where{T} 1.0 end

abstract type ChainSampler{T} end
function sample_chain(::ChainSampler{T})::Vector{T} where {T} end


#abstract type ChainSampler{T} end
#struct MetropolisSampler
