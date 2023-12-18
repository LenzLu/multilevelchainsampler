struct PowerGridProposal <: AbstractProposalGenerator
    net::NetworkEnsemble
    change_size::Int64
    power::Float64 
    coupling::Float64
    damping::Float64
end
PowerGridProposal(net, coupling=1.0, damping=0.1, power=coupling/10.0) = PowerGridProposal(net, 12, power, coupling, damping)

function propose(q::PowerGridProposal)
    G = propose(q.net)
    W = q.power * rand(nv(G)); W .-= mean(W)
    return PowerGrid(G, W, q.coupling, q.damping)
end

function propose(q::PowerGridProposal, g::PowerGrid)
    G = g.grid
    for i=1:q.change_size
        G = propose(q.net, G)
    end
    W = g.power
    #W = q.power * rand(nv(G)); W -= mean(W)
    W = shuffle!(W)
    return PowerGrid(G, W, q.coupling, q.damping)
end