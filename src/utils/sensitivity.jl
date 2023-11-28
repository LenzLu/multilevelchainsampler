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

function sensitivity_inplace(u₀, f!, ∂ᵤf!)
    n = length(u₀)
    v₀ = [u₀ ..., I(n)... ]
    println("v₀ ", length(v₀))

    function F!(dv, v,p,t)
        u = v[1:n]
        G = reshape(v[n+1:n+n^2], n,n)
        J = zeros(eltype(v), n,n)
        f!( dv[1:n], u,p,t )
        ∂ᵤf!( J, u,p,t)
        dv[n+1:n+n^2] .= reshape( G * J, n^2)
    end
    return v₀, F!
end

# TODO: evtl. consider ∂ₜf, ∂ₚf also
