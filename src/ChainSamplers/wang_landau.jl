@Base.kwdef struct EnergyHistogram
    lb = 0.0 ; ub = 1.0
    bins = 100
end

function centers(h::EnergyHistogram)
    x = ([ 1:bins ... ] - .5 ) 
    x .=  (up - lb)/bins; x .+=  lb 
    return x
end

struct WangLandau <: ChainSampler{T}
    histogram::EnergyHistogram
    propose::ProposalGenerator{T}  
end

function sample(algo::WangLandau{T}; x0::T = initialize(algo.propose))
    hist = zeros(algo.histogram.bins)
    
end