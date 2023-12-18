# Metropolis Hastings
@kwdef struct MetropolisHastings <: AbstractChainSampler
    n::Int64 = 100
end


function sample(e::AbstractEnergyFunction, s::MetropolisHastings, g::AbstractProposalGenerator; x0=propose(g))
    E0 = evaluate(e, x0)
    chain = ChainSample([x0], [E0], [0], [(; acceptance=1.0)]) 

    for t=1:s.n
        x = chain.states[end]

        y = propose(g, x)
        q_xy, q_yx = transitions(g, x,y)

        Ex = chain.energies[end]
        Ey = evaluate(e, y)

        A = min(1.0, exp(Ex - Ey) * q_yx/q_xy )

        accept = rand() < A
        if accept
            push!(chain.states, y)
            push!(chain.energies, Ey)
            push!(chain.rejections, 0)
            push!(chain.attributes, (; acceptance=A) )
        else
            chain.rejections[end] += 1
        end    
    end
    return chain
end

# Delayed Acceptance
function sample(e::MultilevelEnergyFunction, s::MetropolisHastings, g::AbstractProposalGenerator; x0=propose(g))
    E0 = zeros(length(e))
    E0[1] = evaluate(e,x0; level=1)
    for i=2:length(e)
        E0[i] = evaluate(e, x0; level=i, cache=E0[i-1]) 
    end
    chain = ChainSample([x0], [E0[end]], [0], 
                        [(; acceptance=1.0, 
                            proxies=Tuple(E0[1:end-1]), 
                            rejections=zeros(Int64, length(e)) )] ) 

    for t=1:s.n
        x = chain.states[end]
        
        y = propose(g, x)
        q_xy, q_yx = transitions(g, x,y)
        
        
        accept=true; A = 1.0;
        proxies = zeros(length(e))
        Ey = nothing; lvl=0 
        while lvl < length(e) ; lvl+=1
            Ex = lvl==length(e) ? chain.energies[end] : chain.attributes[end].proxies[lvl]
            Ey = evaluate(e, y; level=lvl, cache=Ey)
            A = min(1.0, exp(Ex - Ey) * q_yx/q_xy )
            
            if rand() < A
                q_xy = min(1.0, exp(Ey - Ex) * q_xy/q_yx )
                q_yx = A
                proxies[lvl] = Ey
            else
                accept = false
                chain.attributes[end].rejections[lvl] += 1
                break
            end
        end
        if accept
            push!(chain.states, y)
            push!(chain.energies, proxies[end])
            push!(chain.rejections, 0)
            push!(chain.attributes, 
                (; acceptance=A, 
                   proxies = Tuple(proxies[1:end-1]), 
                   rejections=zeros(length(e)) ))
        else
            chain.rejections[end] += 1
        end    
    end
    return chain
end
