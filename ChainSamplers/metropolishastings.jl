include("abstract.jl")

struct MetropolisHastings{T} <: ChainSampler{T}
    energy::Function # Maps T → ℝ
    proposal::ProposalGenerator{T}
    length::Integer
end

function _init_chain(x0, energy::Function, chain_length::Integer)
    T = typeof(x0)

    X = T[x0];        resize!(X, chain_length)
    E = [energy(x0)]; resize!(E, chain_length)
    return X,E
end

function _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)

    A = min(1, exp(Ex - Ey)*q_xy/q_yx)

    # Update
    accept = rand() < A
    return accept ? (y,Ey) : (x,Ex)
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

        X[i+1], E[i+1] = _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)

    end

    return X, E
end

struct DelayedAcceptanceMetropolisHastings{T} <: ChainSampler{T}
    energy::Function # Maps T → ℝ
    surrogates::Vector{Function}
    proposal::ProposalGenerator{T}
    length::Integer
    surrlengths::Vector{<:Integer}

    # TODO: probabilistic subchain length samplers
    #surrlengths::NTuple{Sampleable{Univariate, Discrete}, L}
end


function sample_chain(
    s::DelayedAcceptanceMetropolisHastings;
    x0=initialize(s.proposal))

    L = length(s.surrogates)
    @assert length(s.surrlengths) == L "$L surrogate sampling lengths required"

    # Trivial case
    if L == 0
        subsampler = MetropolisHastings(s.energy, s.proposal, s.length)
        chain = sample_chain(subsampler; x0)
        return chain
    end # from now on L ≧ 1

    # Initialize
    X,E = _init_chain(x0, s.energy, s.length)

    subsampler = DelayedAcceptanceMetropolisHastings(
        s.surrogates[1],  Function[s.surrogates[2:end]...], s.proposal,
        s.surrlengths[1], Integer[s.surrlengths[2:end]...])

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]

        # proposal
        subchain = sample_chain(subsampler, x0=x)
        y = subchain[1][end]
        q_xy = exp(-subchain[2][1])
        q_yx = exp(-subchain[2][end])

        # evaluate energy
        Ey = s.energy(y)

        X[i+1], E[i+1] = _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)
    end

    return X, E
end
