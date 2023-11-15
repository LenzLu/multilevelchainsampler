# include("abstract.jl")
# include("metropolis_hastings.jl")


struct DelayedAcceptanceMetropolisHastings{EnergyType,SurrogateType,T} <:
    ChainSampler{T} where {EnergyType<:Function, SurrogateType<:Function}

    energy::EnergyType
    surrogate::SurrogateType
    proposal::ProposalGenerator{T}
    length::Integer
end


function sample_chain(
    s::DelayedAcceptanceMetropolisHastings;
    x0=initialize(s.proposal))

    X,E = _mh_init_chain(x0, s.energy, s.length)
    F = [s.surrogate(x0)]; resize!(F, s.length) # approximate energies

    subsampler = MetropolisHastings(s.energy, s.proposal, 1)

    for i=1:s.length-1

        # load previous
        x = X[i]; Ex = E[i]; Fx = F[i]

        # proposal
        y = deepcopy(x); propose!(s.proposal, y)
        q_xy,q_yx = transitions(s.proposal, x, y)

        Fy = s.surrogate(y)

        # Promotion step
        g_xy = min(1, exp(Fx - Fy)*q_xy/q_yx)
        if rand() > g_xy
            X[i+1],E[i+1],F[i+1] = x,Ex,Fx # reject early, just copy values
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
