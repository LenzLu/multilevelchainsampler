#Proportional to unnormalized negative logpd, uniform by default
abstract type AbstractEnergyFunction end 
evaluate(::AbstractEnergyFunction, x; args ...  )::Float64 = 0.0

struct SimpleEnergyFunction{F <: Function} <: AbstractEnergyFunction 
    f::F
end 
evaluate(e::SimpleEnergyFunction, x; args...) = e.f(x) # consume arguments


include("multilevel.jl")
include("tracked.jl")


function EnergyFunction(funcs::Vararg{Function}; tracked=false)
    if length(funcs) == 1 
        E = SimpleEnergyFunction(funcs[1]) 
        if tracked
            E = TrackedSimpleEnergyFunction{typeof(E)}(E, 0, 0.0)
        end
    else
        E = ProxiedEnergyFunction([SimpleEnergyFunction.(funcs) ... ])
        if tracked
            E = TrackedMultilevelEnergyFunction(E, zeros(Int,length(E)), zeros(length(E)) )
        end
    end
    return E
end
