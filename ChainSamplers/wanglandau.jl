include("abstract.jl")
include("histutils.jl")

function is_flat(H::Vector, tol=0.95)
    minimum(H) >= tol * mean(H)
end


struct WangLandau{T} <: ChainSampler{T}
    energy::Function # Maps T → ℝ
    proposal::ProposalGenerator{T}
    histogram::EnergyHistogram
    length::Integer
    #entropy_update::Function
end

function sample_chain(
    sampler::WangLandau;
    x0=initialize(sampler.proposal))

    T = typeof(x0)
    X = T[x0]; E = [energy(x0)]
    resize!(X, sampler.length)
    resize!(E, sampler.length)

    # Histogram and entropy
    h = sampler.histogram
    H = zeros(UInt, h.N)
    S = zeros(Float32, h.N)
    f = 1.0

    n = which_bin(h, E[end])
    H[n] += 1; S[n] += f
    N = [n]


    for i=1:sampler.length-1
        x = X[i]; Ex = E[i] # previous

        y = deepcopy(x) # proposal
        propose!(sampler.proposal, y)
        q_xy = probability(sampler.proposal, x, y)
        q_yx = probability(sampler.proposal, y, x)
        Ey = energy(y)

        nx = N[end]
        ny = which_bin(h, Ey)
        A = min(1, exp(S[nx] - S[ny])*q_xy/q_yx)

        accept = rand() < A
        X[i+1] = (accept ?  y :  x)
        E[i+1] = (accept ? Ey : Ex)

        n = accept ? ny : nx
        append!(N, n)
        H[n] += 1; S[n] += f

        if is_flat(H)
            H[:] .= 0
            f *= 0.5  # Refine the f parameter
        end
    end

    m = midpoints(h)
    print("f end ", f)
    plot(m, S, label="entropy")
    display(Plots.current())

    return X, E
end