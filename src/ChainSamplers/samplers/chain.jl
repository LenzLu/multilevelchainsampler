struct ChainSample
    states::Vector         
    energies::Vector
    rejections::Vector
    attributes::Vector
end
Base.length(s::ChainSample) = length(s.states)

function unfold(s::ChainSample)
    states = vcat( map(u -> repeat([u[1]], u[2]), zip(s.states, 1 .+ s.rejections))...  )
    energies = vcat( map(u -> repeat([u[1]], u[2]), zip(s.energies, 1 .+ s.rejections))...  )
    attributes = vcat( map(u -> repeat([u[1]], u[2]), zip(s.attributes, 1 .+ s.rejections))...  )
    return (; states, energies, attributes)
end

Base.getindex(s::ChainSample, attr::Symbol) = getproperty.(s.attributes, attr)

#=
function (+)(s1::ChainSample, s2::ChainSample)
    ChainSample(vcat(s1.states,     s2.states), 
                vcat(s1.energies,   s2.energies), 
                vcat(s1.rejections, s2.rejections))
end
=#

function append!(s1::ChainSample, s2::ChainSample)
    append!(s1.states,     s2.states)
    append!(s1.energies,   s2.energies)
    append!(s1.rejections, s2.rejections)
    append!(s1.attributes, s2.attributes)
end

#=
function getindex(s::ChainSample, i::Int64)
    indices = cumsum(s.repetitions .+ 1) 
    idx = argmax(i .< indices)
    return states[idx], energies[idx]    
end
=#

function save_chain(s::ChainSample, filename::String)
    if !occursin(".", filename) || split(filename, ".")[end] != "jld"
        @warn "Please save as a .jld file instead."
    end

    data = Dict(
        "states" => s.states, 
        "energies" => s.energies,
        "rejections" => s.rejections        
    )
    fields = eltype(s.attributes).parameters[1]
    for field in fields 
        data[string(field)] = s[field]
    end
    data["names"] = [keys(data) ... ]
    save(filename, data)
end
function load_chain(filename::String)
    if !occursin(".", filename) || split(filename, ".")[end] != "jld"
        @warn "Loading a non .jld file."
    end

    data = Dict(
        f => load(filename, f)
        for f in load(filename, "names")
    )
    states = pop!(data, "states")
    energies = pop!(data, "energies")
    rejections = pop!(data, "rejections")
    #attrs = (; (Symbol(k)=>v for (k,v) in data)... )
    attrs = [(; (Symbol(k)=>v[i] for (k,v) in data)...) for i=1:length(states)]
    return ChainSample(states, energies, rejections, attrs)
end


