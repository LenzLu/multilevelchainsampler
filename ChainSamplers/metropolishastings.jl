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
    Pf = [1.0];       resize!(Pf, chain_length) # forward transition probability
    Pb = [1.0];       resize!(Pb, chain_length) # backward transition probabity
    return X,E,Pf,Pb
end

function _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)
    # Acceptance probability
    A = min(1, exp(Ex - Ey)*q_xy/q_yx)
    B = min(1, exp(Ey - Ex)*q_yx/q_xy)

    # Update
    accept = rand() < A
    return accept ? (y,Ey,A,B) : (x,Ex,1-A,1-B)
end

function sample_chain(
    s::MetropolisHastings;
    x0=initialize(s.proposal),
    return_probabilities::Bool = false)

    X,E,Pf,Pb = _init_chain(x0, s.energy, s.length)

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

        X[i+1], E[i+1], Pf[i+1], Pb[i+1] = _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)

    end

    if return_probabilities return X,E,Pf,Pb end

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
    x0=initialize(s.proposal),
    return_probabilities::Bool = false)
    L = length(s.surrogates)

    @assert length(s.surrlengths) == L "$L surrogate sampling lengths required"

    # Trivial case
    if L == 0
        subsampler = MetropolisHastings(s.energy, s.proposal, s.length)
        chain = sample_chain(subsampler; x0, return_probabilities)
        return chain
    end # from now on L ≧ 1

    # Initialize
    X,E,Pf,Pb = _init_chain(x0, s.energy, s.length)

    subsampler = DelayedAcceptanceMetropolisHastings(
        s.surrogates[1],  Function[s.surrogates[2:end]...], s.proposal,
        s.surrlengths[1], Integer[s.surrlengths[2:end]...])

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]

        # proposal
        subchain = sample_chain(subsampler, x0=x, return_probabilities=true)
        y = subchain[1][end]
        a,b = subchain[3], subchain[4]
        q_xy = prod( a./b )
        q_yx = prod( b./a )

        # evaluate energy
        Ey = s.energy(y)

        X[i+1], E[i+1], Pf[i+1], Pb[i+1] = _acceptance_step(x,y,Ex,Ey,q_xy,q_yx)
    end

    if return_probabilities return X,E,Pf,Pb end

    return X, E
end
