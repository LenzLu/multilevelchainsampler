include("abstract.jl")

# Connected erdoes renyi
struct ErdoesRenyiSampler <: NetworkSampler
    N::Integer
    p::Real # 0 <= p <= 1
end

function initialize(S::ErdoesRenyiSampler)
    G = erdos_renyi(S.N, S.p)
    while !is_connected(G)
        G = erdos_renyi(S.N, S.p)
    end
    return G
end

function propose!(S::ErdoesRenyiSampler, G::Graph)
    u,v = choose_node_pair(G)
    if has_edge(G, u, v)
        rem_edge!(G, u, v)
    else
        if rand() < S.p
            add_edge!(G, u,v)
        end
    end
end
