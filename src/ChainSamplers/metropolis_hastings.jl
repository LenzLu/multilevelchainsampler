# include("abstract.jl")

struct MetropolisHastings{EnergyType,T} <: ChainSampler{T} where {EnergyType <: Function}
    energy::EnergyType
    proposal::ProposalGenerator{T}
    length::Int64
end


function _mh_init_chain(x0, energy::Function, chain_length::Integer)
    T = typeof(x0)
    X = T[x0];        resize!(X, chain_length)
    E = [energy(x0)]; resize!(E, chain_length)
    return X,E
end

function sample_chain(
    s::MetropolisHastings;
    x0=initialize(s.proposal))

    X,E = _mh_init_chain(x0, s.energy, s.length)

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]

        # proposal
        y = deepcopy(x); propose!(s.proposal, y)
        q_xy,q_yx = transitions(s.proposal, x, y)

        # evaluate energy
        Ey = s.energy(y)

        # update
        A = min(1, exp(Ex - Ey)*q_xy/q_yx)
        accept = rand() < A
        X[i+1],E[i+1] = accept ? (y,Ey) : (x,Ex)
    end

    return X, E
end
