# TODO: templates

struct PowerGrid
    grid
    power
    coupling 
    damping
    syncstate
end

function PowerGrid(grid, power, coupling=1.0, damping=0.1) 
    syncstate = synchronous_state(grid, power, coupling, damping)
    g =PowerGrid( grid, power, coupling, damping, syncstate )
    

    if !stable(g, randn(2*nv(grid))) 
        @warn "Synchronous state not found!" end
    return g
end
