using Plots

include("dynamics.jl")


# Distribution Histograms
x = sample_perturbation(5000)
stephist(u₀[1] .+ x[:,1], alpha=.25, normalize=true, width=1, color=:black, label="Initial")

blend(c1, c2, t) = c1 .* (1-t) .+ c2 .* t

timesteps = [0.5, 0.33, 0.25, 0.1, 0.05, 0.01]
for (i,dt) = enumerate(timesteps)
    q = quantity_of_interest(x, dt)
    n = length(timesteps); color = RGB(i/n, 1-i/n, 1-i/n)
    stephist!(q, alpha=.2, normalized=true, label="dt=$dt", width=4; color)
end

title!("Distribution mapping")
imgdir = "$(@__DIR__)/imgs"
savefig("$imgdir/distributions.png")
display(Plots.current())

# Statistics
timesteps = 0.1 .^ [0:.1:3 ... ]
plot()
xlims!(extrema(timesteps)... )
nsamples = [2,10]; K = 10
for (i,N) = enumerate(nsamples)
    m = zeros(length(timesteps), K)
    s = zeros(length(timesteps), K)
    for k=1:K
        for (j,dt) in enumerate(timesteps)
            x = sample_perturbation(N)
            q = quantity_of_interest(x, dt)
            m[j,k] = mean(q)
            s[j,k] = std(q)
        end
    end
    n = length(nsamples)
    cm = RGB(blend([1,0,0], [.5,0,0], i/n)...)
    cs = RGB(blend([0,0,1], [0,0,.5], i/n)...)
    errorline!(timesteps, m, color=cm, label="mean N=$N")
    errorline!(timesteps, s[:,:], color=cm, label="std N=$N")
end
xlabel!("Time step Δt")
title!("Statistics")
#xaxis!(:log); yaxis!(:log)
savefig("$imgdir/statistics.png")
display(Plots.current())


# Solver time
timesteps = 0.1 .^ [0:.25:4 ... ]
plot(); xlims!(extrema(timesteps)...); t = []
for N = [5, 10, 50, 100]
    x = sample_perturbation(N); t = []
    for dt in timesteps
        τ = @elapsed quantity_of_interest(x, dt)
        append!(t, τ)
    end
    plot!(timesteps, t, label="N = $N")
end
xaxis!(:log); yaxis!(:log)
xlabel!("Time step Δt")
ylabel!("Computation time")
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
println("Error rate O(hᵞ), γ = $γ ")
plot!(timesteps, C .* timesteps.^γ, label="O(hᵞ), γ = $(round(γ,digits=2)) ")

xaxis!(:log); yaxis!(:log)
xlabel!("Time step Δt")
ylabel!("Relative error")
title!("Convergence rate")
savefig("$imgdir/convergence.png")
display(Plots.current())
