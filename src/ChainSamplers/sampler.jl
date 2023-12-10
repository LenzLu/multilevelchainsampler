
## Proposal generator

abstract type ProposalGenerator{StateType} end
function initialize(::ProposalGenerator{T})::T where{T} zero(T) end
function propose!(x::ProposalGenerator{T}, ::T) where{T} rand!(x) end

# Default transition(s)
function transition(::ProposalGenerator{T}, ::T, ::T) where{T} return 1.0 end
function transitions(g::ProposalGenerator{T}, x::T, y::T) where{T} 
    return (transition(g,x,y), transition(g,y,x)) end


struct ChainSample{X}
    states::Vector{X}
    energies::Vector{Float64}
    repetitions::Vector{Int32}
end

abstract type AbstractChainSampler{X} end
#function sample(s::AbstractChainSamper; x0 = initial)  end
