abstract type MultilevelEnergyFunction <: AbstractEnergyFunction end

struct ProxiedEnergyFunction <: MultilevelEnergyFunction
    proxies::Vector{AbstractEnergyFunction}
end
Base.length(e::ProxiedEnergyFunction) = length(e.proxies)
function evaluate(e::ProxiedEnergyFunction, x; level=length(e), cache=nothing, args...)
    return evaluate(e.proxies[level], x, args...)
end


struct SamplingBasedEneryFunction <: MultilevelEnergyFunction
    xsamples::Vector
    nsamples::Vector{Int64}
    func::Function
end
Base.length(e::SamplingBasedEneryFunction) = length(e.nsamples)
function evaluate(e::SamplingBasedEneryFunction, x; level=length(e), cache=nothing, args...)
    if level==1 || isnothing(cache) 
        E = e.func.( x, e.xsamples[1:e.nsamples[level]])
        return mean( E )
    else
        E = e.func.( x, e.xsamples[e.nsamples[level-1]+1:e.nsamples[level]])
        return cache * (e.nsamples[level-1]/e.nsamples[level]) + sum(E)/e.nsamples[level]
    end
end


