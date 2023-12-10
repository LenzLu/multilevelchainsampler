struct PowerGrid
    grid::Graph
    power::Vector 
    coupling::Vector
    damping::Vector
    syncstate::Vector
end


function PowerGrid(grid::Graph, power::Vector, damping::AbstractFloat, coupling::AbstractFloat) 
    dtype = Float64
    power = convert(Vector{dtype}, power)
    damping = convert(Vector{dtype}, repeat([damping], nv(grid)))
    coupling = convert(Vector{dtype}, repeat([coupling], ne(grid)))
    syncstate = synchronous_state(grid, power, coupling, damping)
    PowerGrid( grid, power, coupling, damping, syncstate )
end