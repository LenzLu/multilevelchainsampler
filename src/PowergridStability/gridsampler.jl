# include("abstract.jl")
# include("../NetworkEnsembles/all.jl")
using Random

## Chain sampling

struct PowerGridEnsemble <: ProposalGenerator{PowerGrid}
    networksampler::NetworkSampler
    change_size::Integer
    K::Real
    α::Real


    #TODO: probabilistic parameters
    #K::Sampleable{Univariate, Continuous}
end


function initialize(e::PowerGridEnsemble)
    G = initialize(e.networksampler)

    W = rand(nv(G)); W = W .- mean(W)
    # n = nv(G) ÷ 2;
    # W = ones(2n); W[1:n] .= -1


    grid = PowerGrid( G,W, e.K,e.α )
    return grid
end

function propose!(ensemble::PowerGridEnsemble, grid::PowerGrid)
    G = grid.G
    for i=1:ensemble.change_size
        propose!(ensemble.networksampler, G)
    end
    W = rand(nv(grid.G)); W = W .- mean(W); grid.W[:] = W

    # n = nv(G) ÷ 2;
    #W = ones(2n); W[1:n] .= -1
    #shuffle!(grid.W)

    return grid
end
