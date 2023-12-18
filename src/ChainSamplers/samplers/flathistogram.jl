struct WangLandau <: AbstractChainSampler
    levels::Vector{Float64} # energy histogram edges
    n::Int64
    entropy::Vector{Float64} 
end
function WangLandau(lb::Number, ub::Number, n::Int64=100; nbins::Integer=100)
    return WangLandau(lb:(ub - lb)/nbins:ub, n, zeros(nbins))
end

function sample(e::AbstractEnergyFunction, s::WangLandau, g::AbstractProposalGenerator; x0=propose(g), learning_rate=1.0)
    E0 = evaluate(e, x0)
    bin_x = searchsortedfirst(s.levels, E0)
    bin_x = max(1, min(length(s.levels)-1, bin_x))

    chain = ChainSample([x0], [E0], [0], [(; acceptance=1.0, bin=bin_x)]) 

    for t=1:s.n
        x = chain.states[end]
        y = propose(g, x)

        Ey = evaluate(e, y)
        bin_y = searchsortedfirst(s.levels, Ey)
        bin_y = max(1, min(length(s.levels)-1, bin_y))

        Sx = s.entropy[bin_x]
        Sy = s.entropy[bin_y]

        q_xy, q_yx = transitions(g, x,y)
        A = min(1.0, exp(Sx - Sy) * q_yx/q_xy )

        accept = rand() < A
        if accept
            push!(chain.states, y)
            push!(chain.energies, Ey)
            push!(chain.rejections, 0)
            push!(chain.attributes, (; acceptance=A, bin=bin_y))
            bin_x = bin_y
        else
            chain.rejections[end] += 1
        end    
        s.entropy[bin_x] += learning_rate / t
    end

    return chain
end    

#=struct WangLandau <: AbstractChainSampler
    levels::Vector{Float64} # energy histogram edges
    entropy::Vector{Float64} 
    flatness::Float64 
    temperature::Vector{Float64}
    stop::Float64
    step::Function 
end
function WangLandau(lb::Number, ub::Number,
    nbins::Integer=100; flatness = 0.8, 
    temperature = [exp(1.0)], stop = 1e-6,
    step = ((T,t) -> 1/t))
    levels = lb:(ub - lb)/nbins:ub
    entropy = ones(nbins)
    return WangLandau(levels, entropy, flatness, temperature, stop, step)
end

function is_flat(histogram::Vector{Int64}, use_nonzero=false, flatness=0.6)
    h = histogram ./ sum(histogram)
    if (use_nonzero) h = h[(h .> 0.0)] end
    return ( minimum(h) > flatness * mean(h) ) 
end 

function sample(e::AbstractEnergyFunction, s::WangLandau, g::AbstractProposalGenerator; x0=propose(g), notify_interval=1000) 
    E0 = evaluate(e, x0)
    chain = ChainSample([x0], [E0], [0]) 

    histogram = zeros(Int64, length(s.entropy))
    bin_x = searchsortedfirst(s.levels, E0)
    bin_x = max(1, min(length(s.levels)-1, bin_x))

    t = 0
    while s.temperature[end] > s.stop
        x = chain.states[end]
        y = propose(g, x)

        Ey = evaluate(e, y)
        bin_y = searchsortedfirst(s.levels, Ey)
        bin_y = max(1, min(length(s.levels)-1, bin_y))

        Sx = s.entropy[bin_x]
        Sy = s.entropy[bin_y]

        q_xy, q_yx = transitions(g, x,y)
        A = min(1.0, exp(Sx - Sy) * q_yx/q_xy )

        accept = rand() < A
        if accept
            push!(chain.states, y)
            push!(chain.energies, Ey)
            push!(chain.rejections, 0)
            bin_x = bin_y
        else
            chain.rejections[end] += 1
        end    
        histogram[bin_x] += 1
        s.entropy[bin_x] += s.temperature[end]

        if  is_flat(histogram, t>1000, s.flatness)
            histogram .= 0
            histogram[bin_x] += 1
            T = s.step(s.temperature[end], t)
            push!(s.temperature, T)
        end

        t += 1
        #=if t % notify_interval == 0 
            println("Iteration ", t, " temperature ", s.temperature)
            println("Histogram ", histogram)
        end=#
    end
    return chain
end
=#