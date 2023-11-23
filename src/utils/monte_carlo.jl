using Statistics

function estimator(sampler::Function, quantity::Function, nsamples::Int64=128;
    h=0.1, antithetic=true)

    x = sampler(nsamples)
    if antithetic x = vcat( x, .-x) end
    Q = quantity(x, h)
    return mean(Q), std(Q)
end

function multilevel_estimator(
    sampler::Function,
    quantity::Function,
    nsamples::Vector{Int64} = [128,64,32,16,8];
    h0 = 0.1, m=0.5 )

    # Coursest level
    x = sampler(nsamples[1])
    Q₀ = quantity(x, h0)
    Y  = mean(Q₀)
    s² = std(Q₀)^2

    # Iterate over higher levels
    for l in 2:length(nsamples)
        x = sampler(nsamples[l])
        hₗ, hₗ₋₁ = @. h0*m^( [l, l-1] -1)
        Qₗ   = quantity(x, hₗ  )
        Qₗ₋₁ = quantity(x, hₗ₋₁)
        Y  += mean(Qₗ .- Qₗ₋₁)
        s² += std( Qₗ .- Qₗ₋₁)^2
    end

    return Y, sqrt(s²)
end

function multilevel_estimator(
    sampler::Function,
    quantities::FunctionIterable;
    nsamples::Vector{Int64} = 4 .^ [length(quantities):-1:1 ... ] )
    @assert length(nsamples) == length(quantities) "sampling levels incoherent, give $(length(quantities)) integers"

    # Coursest level
    x = sampler(nsamples[1])
    Q₀ = quantities[1](x)
    Y  = mean(Q₀)
    s² = std(Q₀)^2

    # Iterate over higher levels
    for l in 2:length(nsamples)
        x = sampler(nsamples[l])
        Qₗ   = quantities[l-1](x)
        Qₗ₋₁ = quantities[l](x)
        Y  += mean(Qₗ .- Qₗ₋₁)
        s² += std( Qₗ .- Qₗ₋₁)^2
    end

    return Y, sqrt(s²)
end

#=
function adaptive_multilevel_estimator(
    sampler::Function, quantity::Function;
    h0 = 1.0, m = 0.5,
    tolerance=1e-3,
    n0=10, α=1.0, r=1.0)

    N = repeat([n0], 2)
    Y = [[], []]
    costs = []

    cost_estimate(level) = (@elapsed quantity(sampler(n0), h0*m^(level-1)))

    finished = false
    while !finished
        L = length(N)

        #estimate costs with single random sample
        costs = [costs..., ( cost_estimate(l) for l=length(costs)+1:L ) ... ]

        # update samples
        x  = sampler(N[1] - length(Y[1]))
        Q₀ = quantity(x, h0)
        append!(Y[1], Q₀)

        for l=2:L
            x = sampler(N[l] - length(Y[l]))
            hₗ, hₗ₋₁ = @. h0*m^( [l, l-1] -1)
            Qₗ   = quantity(x, hₗ  )
            Qₗ₋₁ = quantity(x, hₗ₋₁)
            append!(Y[l], Qₗ - Qₗ₋₁)
        end

        # update estimators
        μ = mean.(Y); σ = std.(Y)

        N = @. max(N, convert(Int64,ceil(σ/sqrt(costs))) )
        bound = ((r*m^α - 1.0)/sqrt(2.0))  * tolerance
        if abs(μ[L]) > abs(bound)
            append!(N, [n0])
            append!(Y, [[]])
            print("L = ", L+1)
        else
            err =  sum( σ.^2 ./ N ) - tolerance^2/2
            println("error ", err, " N ", N, " C ", costs, " b ", bound)
            if sum( σ.^2 ./ N ) <= tolerance^2/2
                finished = true
            end
        end
    end

    return sum(μ) , sqrt( sum(σ.^2) )
end
=#
