struct PowerGrid
    grid
    power::Vector
    coupling
    damping
    syncstate::Vector
end


function swing_dynamics!(dₜx, x, p, t)
    G,W,K,α = p.grid, p.power, p.coupling, p.damping 
    N = nv(G)
    
    dₜx[N+1:2N] .= W .- α .* x[N+1:2N]
    for (k,e) in enumerate(edges(G))
        i,j = e.dst, e.src
        x_ij = x[i] - x[j]

        dₜx[N+i] -= K * sin(x_ij)
        dₜx[N+j] += K * sin(x_ij)
    end
    dₜx[1:N] .= x[N+1:2N];
end
