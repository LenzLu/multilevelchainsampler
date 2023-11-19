using LinearAlgebra
using OrdinaryDiffEq

function sensitivity(u₀, f, ∂ᵤf)
    n = length(u₀)
    println("n = $n")
    v₀ = [u₀ ..., I(n)... ]
    #v₀ = [u₀ ... , I(n)... , -f(u₀, p, t₀)... ]

    function F(v,p,t)
        u = v[1:n]
        G = reshape(v[n+1:n+n^2], n,n)

        dₜu = f(u,p,t)
        dₜG = G*reshape(∂ᵤf(u,p,t),n,n)

        dₜv = [dₜu ... , dₜG ... ]
        return dₜv
    end

    return v₀, F
end

# Caveat: Inplace evaluation tricky since ∂ᵤf! would change state

# TODO: evtl. consider ∂ₜf, ∂ₚf also
