include("dynamics.jl")

using Plots
using Revise
using MultilevelChainSampler




# Solver time
timesteps = 0.1 .^ [0:.25:4 ... ]
plot(); xlims!(extrema(timesteps)...); t = []
for N = [10, 50, 100, 150]
    x = sample_perturbation(N); t = []
    for dt in timesteps
        τ = @elapsed quantity_of_interest(x, dt)
        append!(t, τ)
    end
    plot!(timesteps, t, label="N = $N")
end
xlabel!("Time step Δt")
ylabel!("Computation time")
title!("Computation cost O(hᵅ), α = $α ")
imgdir = "$(@__DIR__)/imgs"xaxis!(:log); yaxis!(:log)
savefig("$imgdir/solvertime.png")
display(Plots.current())

x = -log.(timesteps)
y = -log.(t)
A = hcat( ones(length(x)), x  )
_,α = inv(A'A)*A'y
xaxis!(:log); yaxis!(:log)
println("Computation cost O(hᵅ), α = $α ")



# Solver accuracy
timesteps = 0.1 .^ [0:.25:4 ... ]
plot(); xlims!(extrema(timesteps)...); t = []
q_ref = quantity_of_interest(zeros(1,2), 1e-5;
    solver=Rodas42(), adaptive=true, reltol=1e-6, abstol=1e-6)[1]
err = []
for dt in timesteps
    q = quantity_of_interest(zeros(1, 2), dt)[1]
    e = abs(q - q_ref)/abs(q_ref)
    append!(err, e)
end
plot!(timesteps, err, label="measured")

x = log.(timesteps[2:end])
y = log.(err[2:end])
A = hcat( ones(length(x)), x  )
C,γ = inv(A'A)*A'y
println("Linear fit")
plot!(timesteps, C .* timesteps.^γ, label="O(hᵞ), γ = $(round(γ,digits=2)) ")

xaxis!(:log); yaxis!(:log)
xlabel!("Time step Δt")
ylabel!("Relative error")
title!("Consistency error O(hᵞ), γ = $γ ")
savefig("$imgdir/consistency.png")
display(Plots.current())
