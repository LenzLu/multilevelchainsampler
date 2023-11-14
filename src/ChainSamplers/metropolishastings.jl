include("abstract.jl")

struct MetropolisHastings{EnergyType,T} <:
    ChainSampler{T} where {EnergyType <: Function, n<:Integer}

    energy::EnergyType                           # Maps T → ℝ
    proposal::ProposalGenerator{T}
    length::Int64
end

struct DelayedAcceptanceMetropolisHastings{EnergyType,SurrogateType,T} <:
    ChainSampler{T} where {EnergyType<:Function, SurrogateType<:Function}

    energy::EnergyType
    surrogate::SurrogateType
    proposal::ProposalGenerator{T}
    length::Integer
end

#=const=# FunctionIterable = Union{AbstractVector{<:Function}, Tuple{Vararg{Function}}}
struct MultilevelMetropolisHastings{EnergyType,SurrogateType,T} <:
    ChainSampler{T} where {EnergyType <: Function, SurrogateType <: FunctionIterable}

    energy::EnergyType
    surrogates::SurrogateType
    proposal::ProposalGenerator{T}
    length::Integer
    surrlengths::Vector{Int64}

    # TODO: probabilistic subchain length samplers
    #surrlengths::NTuple{Sampleable{Univariate, Discrete}, L}

#    @assert length(surrlengths) >= length(surrogates)
#            "At least $L surrogate sampling lengths required"
end



function _init_chain(x0, energy::Function, chain_length::Integer)
    T = typeof(x0)
    X = T[x0];        resize!(X, chain_length)
    E = [energy(x0)]; resize!(E, chain_length)
    return X,E
end

function sample_chain(
    s::MetropolisHastings;
    x0=initialize(s.proposal))

    X,E = _init_chain(x0, s.energy, s.length)

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]

        # proposal
        y = deepcopy(x)
        propose!(s.proposal, y)
        q_xy = probability(s.proposal, x, y)
        q_yx = probability(s.proposal, y, x)

        # evaluate energy
        Ey = s.energy(y)

        # update
        A = min(1, exp(Ex - Ey)*q_xy/q_yx)
        accept = rand() < A
        X[i+1],E[i+1] = accept ? (y,Ey) : (x,Ex)
    end

    return X, E
end


function sample_chain(
    s::DelayedAcceptanceMetropolisHastings;
    x0=initialize(s.proposal))

    X,E = _init_chain(x0, s.energy, s.length)
    F = [s.surrogate(x0)]; resize!(F, s.length)

    subsampler = MetropolisHastings(s.energy, s.proposal, 1)

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]; Fx = F[i]

        # proposal
        y = deepcopy(x)
        propose!(s.proposal, y)
        q_xy = probability(s.proposal, x, y)
        q_yx = probability(s.proposal, y, x)

        Fy = s.surrogate(y)

        # Propomotion step
        g_xy = min(1, exp(Fx - Fy)*q_xy/q_yx)
        if rand() > g_xy
            X[i+1],E[i+1],F[i+1] = x,Ex,Fx # reject early
        else
            g_yx = min(1, exp(Fy - Fx)*q_yx/q_xy)

            # evaluate energy
            Ey = s.energy(y)

            # update
            A = min(1, exp(Ex - Ey) * (q_xy * g_yx) / (q_yx * g_xy) )
            accept = rand() < A
            X[i+1],E[i+1],F[i+1] = accept ? (y,Ey,Fy) : (x,Ex,Fx)
        end

    end

    return X, E
end


function sample_chain(
    s::MultilevelMetropolisHastings;
    x0=initialize(s.proposal),
    return_probabilities::Bool = false)

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

    X,E = _init_chain(x0, s.energy, s.length)

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
        A = min(1, exp(Ex - Ey)*exp(Fy - Fx))
        accept = rand() < A
        X[i+1],E[i+1] = accept ? (y,Ey) : (x,Ex)
    end

    return X, E
end
