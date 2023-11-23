# include("abstract.jl")
function choose_node_pair(G::Graph)
    u = rand(1:nv(G))
    v = rand(1:nv(G)-1)
    if (v >= u) v += 1 end
    return u,v
end

# Connected erdoes renyi
struct ErdoesRenyiSampler <: NetworkSampler
    N::Integer
    p::Real # 1/N < p <= 1
end

function initialize(S::ErdoesRenyiSampler)
    G = erdos_renyi(S.N, S.p)
    while !is_connected(G) #TODO: check if p ≧ log(N)/N
        G = erdos_renyi(S.N, S.p)
    end
    return G
end

function propose!(S::ErdoesRenyiSampler, G::Graph)
    @assert is_connected(G) "only for connected graphs"

    finished = false
    while !finished
        u,v = choose_node_pair(G)
        if has_edge(G, u, v)
            if rem_edge!(G, u, v)
                if is_connected(G)
                    finished = true
                else
                    add_edge!(G,u,v)
                end
            end
        else
            if rand() < S.p
                if add_edge!(G, u,v)
                    finished  = true
                end
            end
        end
    end
end
