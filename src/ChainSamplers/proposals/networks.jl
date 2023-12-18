abstract type NetworkEnsemble <: AbstractProposalGenerator end

## 
struct ErdoesRenyiEnsemble <: NetworkEnsemble
    N::Int64
    p::Float64
    seed::Int64 
end
ErdoesRenyiEnsemble(N,p=log(N)/N * 1.2) = ErdoesRenyiEnsemble(N,p,rand(Int))


function propose(R::ErdoesRenyiEnsemble)          #TODO: check if p ≧ log(N)/N
    G = erdos_renyi(R.N, R.p, seed=R.seed)
    while !is_connected(G)                  
        G = erdos_renyi(R.N, R.p)
    end
    return G
end

function propose(R::ErdoesRenyiEnsemble, G::Graph)
    @assert is_connected(G) "only for connected graphs"
    H = deepcopy(G)

    finished = false
    while !finished
        u,v = samplepair(R.N)
        if has_edge(H, u, v)
            if rem_edge!(H, u, v)
                if is_connected(H)
                    finished = true
                else
                    add_edge!(H,u,v)
                end
            end
        else
            if rand() < R.p
                if add_edge!(H, u,v)
                    finished  = true
                end
            end
        end
    end

    return H
end

struct WattsStrogatzEnsemble <: NetworkEnsemble
    N::Int64
    k::Integer
    p::Float64   # 0 <= p <= 1
    seed::Int64
end
WattsStrogatzEnsemble(N,k=2,p=log(N)/N * 1.2) = WattsStrogatzEnsemble(rand(Int))


function propose(R::WattsStrogatzEnsemble)
    G = watts_strogatz(R.N, R.k, R.p, seed=R.seed)
    while !is_connected(G)
        G = watts_strogatz(R.N, R.k, R.p)
    end
    return G
end

function propose!(R::WattsStrogatzEnsemble, G::Graph)
    H = deepcopy(G)
    E = [edges(H)...]
    finished = false
    while !finished

        # sample edge
        e = E[rand(1:ne(H)) ]
        u,v,w = e.src, e.dst, 0

        if rand() < R.p
            w = mod(u +  (-1)^rand(0:1) * rand(1:R.k), nv(H))
            if (w == 0) w = nv(G) end
        else
            w = rand(1:nv(H)-1)
            if (w >= u) w += 1 end
        end

        if !has_edge(u, w)
            rem_edge!(H, u,v)
            add_edge!(H, u,w)
            finished = true
        end
    end

    return H
end
