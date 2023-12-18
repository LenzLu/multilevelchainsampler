function r_uni(g::PowerGrid)
    x = g.syncstate
    r = 0.0
    for e in edges(g.grid)
        i,j = e.dst, e.src
        r += g.coupling * cos(x[i] - x[j])
    end
    return r / ne(g.grid)
end

function aspl(g::PowerGrid)
    D = johnson_shortest_paths(g.grid)
    return mean(D.dists)
end