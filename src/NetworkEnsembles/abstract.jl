# include("../ChainSamplers/abstract.jl")

using Graphs

const NetworkSampler = ProposalGenerator{Graph}

function choose_node_pair(G::Graph)
    u = rand(1:nv(G))
    v = rand(1:nv(G)-1)
    if (v >= u) v += 1 end
    return u,v
end
