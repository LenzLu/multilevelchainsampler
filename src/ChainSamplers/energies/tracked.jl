#function TrackedEnergyFunction()::AbstractEnergyFunction end

mutable struct TrackedSimpleEnergyFunction{F <: AbstractEnergyFunction} <: AbstractEnergyFunction
    f::F
    evaluations::Int64
    avgtime::Float64
end
TrackedEnergyFunction(f::AbstractEnergyFunction) = TrackedEnergyFunction{typeof(f)}(f, 0, 0.0)

struct TrackedMultilevelEnergyFunction{F <: MultilevelEnergyFunction} <: MultilevelEnergyFunction
    f::F
    evaluations::Vector{Int64}
    avgtime::Vector{Float64}
end
Base.length(e::TrackedMultilevelEnergyFunction) = length(e.f)
TrackedEnergyFunction(f::MultilevelEnergyFunction) = TrackedMultilevelEnergyFunction{typeof(f)}(f, zeros(Int, length(f)), zeros(length(f)))


function evaluate(e::TrackedSimpleEnergyFunction, x; args...) 
    Δt = @elapsed Ex = evaluate(e.f, x; args...)
    e.avgtime = ( e.evaluations * e.avgtime + Δt ) / (e.evaluations+1) 
    e.evaluations += 1
    return Ex
end

function evaluate(e::TrackedMultilevelEnergyFunction, x; args...)
    if !hasproperty(args, :level) 
        args = (; level=length(e), args...) 
    end
    Δt = @elapsed Ex = evaluate(e.f, x; args...)
    e.avgtime[args.level] = ( e.evaluations[args.level] * e.avgtime[args.level] + Δt ) / (e.evaluations[args.level]+1) 
    e.evaluations[args.level] += 1
    return Ex
end
