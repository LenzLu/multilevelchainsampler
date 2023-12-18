## Proposal generator
abstract type AbstractProposalGenerator end

# Random walk by default
function propose(g::AbstractProposalGenerator)    return 0.0 end
function propose(g::AbstractProposalGenerator, x) return x + randn() end

# Symetric transition(s) by default
transition(g::AbstractProposalGenerator, x,y)::Float64 = 1.0
transitions(g::AbstractProposalGenerator, x,y)::NTuple{2, Float64} =
    (transition(g,x,y), transition(g,y,x)) 


include("walk.jl")
include("networks.jl")