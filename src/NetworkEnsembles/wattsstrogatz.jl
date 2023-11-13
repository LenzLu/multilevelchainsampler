include("abstract.jl")

# Connected erdoes renyi
struct WattsStrogatzSampler <: NetworkSampler
    N::Integer
    k::Integer
    p::Real # 0 <= p <= 1
end

function initialize(S::WattsStrogatzSampler)
    G = watts_strogatz(S.N, S.k, S.p)
    while !is_connected(G)
        G = watts_strogatz(S.N, S.k, S.p)
    end
    return G
end

function propose!(S::WattsStrogatzSampler, G::Graph)
    E = [edges(G)...]
    finished = false
    while !finished

        # sample edge
        e = E[rand(1:ne(G)) ]
        u,v = e.src, e.dst

        if rand() < S.p
            w = mod(u +  (-1)^rand(0:1) * rand(1:S.k), nv(G))
            if (w == 0) w = nv(G) end
        else
            w = rand(1:nv(G)-1)
            if (w >= u) w += 1 end
        end

        if !has_edge(u, w)
            rem_edge!(G, u,v)
            add_edge!(G, u,w)
            finished = true
        end
    end
end
