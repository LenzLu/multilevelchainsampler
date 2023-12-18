function sync_rhs(x,p)
    N = length(x)÷2
    dₜx = similar(x)
    swing_dynamics!(dₜx, x, p, 0.0)
    dₜx[1:N] .-= mean(dₜx[1:N])
    return dₜx[2:2N]
end

function synchronous_state(grid, power, coupling, damping)
    p = (; grid,power,coupling,damping)
    N = nv(grid)

    #= 
    # Linear stability initialization
    u0 = zeros(2N)    
    L = Matrix(laplacian_matrix(grid))    
    U,S,V = svd(L)
    z = abs.(S) .< 1e-9; nz = .! z
    S[z] .= 0; S[nz] .= 1 ./ S[nz]
    ϕ0 = V * diagm(S) * U' * power 
    #ϕ0 = mod.(ϕ0 .+ π, 2π) .-π
    println("ϕ0")
    u0[1:N] .= ϕ0 ./ coupling
    =#

    nl = NonlinearProblem(sync_rhs, zeros(2N), p)
    sync = solve(nl, NewtonRaphson()).u
    #sync[1:N] .-= mean(sync[1:N])
    return sync
end 

function PowerGrid(grid, power::Vector, coupling=1.0, damping=0.1) 
    syncstate = synchronous_state(grid, power, coupling, damping)
    g = PowerGrid( grid, power, coupling, damping, syncstate )

    if !stable(g, zeros(2*nv(grid))) 
        @warn "Synchronous state not found!" end
    return g
end
