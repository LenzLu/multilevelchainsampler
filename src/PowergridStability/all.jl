using Graphs

struct PowerGrid
    G::Graph
    W::Vector
    K::Real
    α::Real
end

include("gridsampler.jl")
include("dynamics.jl")
#include("display_gplots.jl")
include("perturbations.jl")
include("stability.jl")
include("stability_sens_proxy.jl")
