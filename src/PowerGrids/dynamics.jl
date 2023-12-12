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

function sync_rhs(x,p)
    N = length(x)÷2
    dₜx = similar(x)
    swing_dynamics!(dₜx, x, p, 0.0)
    dₜx[1:N] .-= mean(dₜx[1:N])
    return dₜx[2:2N]
end


function synchronous_state(grid, power, coupling=1.0, damping=0.1)
    p = (; grid,power,coupling,damping)
    N = nv(grid)
    nl = NonlinearProblem(sync_rhs, zeros(2N), p)
    sync = solve(nl, NewtonRaphson()).u
    sync[1:N] .-= mean(sync[1:N])
    return sync
end 
