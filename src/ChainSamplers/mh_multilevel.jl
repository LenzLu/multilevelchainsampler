
struct MultilevelMetropolisHastings{EnergyType,SurrogateType,T} <:
    ChainSampler{T} where {EnergyType<:Function, SurrogateType<:FunctionIterable}

    energy::EnergyType
    surrogates::SurrogateType
    proposal::ProposalGenerator{T}
    length::Integer
    surrlengths::Vector{Int64}
end


function sample_chain(
    s::MultilevelMetropolisHastings;
    x0=initialize(s.proposal))

    # Trivial case
    L = length(s.surrogates)
    if L == 0
        subsampler = MetropolisHastings(s.energy, s.proposal, s.length)
        chain = sample_chain(subsampler; x0)
        return chain
    end # from now on L ≧ 1

    # Recursion
    subsampler = MultilevelMetropolisHastings(
        s.surrogates[1],  s.surrogates[2:end], s.proposal,
        s.surrlengths[1], s.surrlengths[2:end] )

    X,E = _mh_init_chain(x0, s.energy, s.length)

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]

        # proposal
        subchain = sample_chain(subsampler, x0=x)
        y = subchain[1][end]
        Fx,Fy = subchain[2][1], subchain[2][end]

        # evaluate energy
        Ey = s.energy(y)

        # update
        A = min(1, exp( (Ex - Ey) + (Fy - Fx) ) )
        accept = rand() < A
        X[i+1],E[i+1] = accept ? (y,Ey) : (x,Ex)
    end

    return X, E
end
