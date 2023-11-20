include("dynamics.jl")

using Plots

v0 = deviation(zeros(2); solver=RK4(), dt=1e-3)
q0 = v0[1, end]
plot(v0[1,:], v0[2,:], label="unperturbed", color="black")

δu = sample_perturbation(100)
n = size(δu,1)
V = [deviation(δu[i,:]; solver=RK4(), dt=1e-3) for i=1:n]
Q = [v[1,end] for v in V]

for (i,v) in enumerate(V)
    plot!(v[1,:], v[2,:], label="", color=RGBA(.2,.2,.3,.1))
end
m = u₀ .+ mean(δu, dims=1)[1,:]
s = std(δu, dims=1)[1,:]
plot!(m[1] .+ s[1] .* [1,1,-1,-1,1], m[2] .+ s[2] .* [1,-1,-1,1,1], color="black", linestyle=:dot, label="initial perturbation")

m = [ mean([v[i,end] for v in V]) for i=[1,2] ]
s = [ std([v[1,end] for v in V]) for i=[1,2] ]
plot!(m[1] .+ s[1] .* [1,1,-1,-1,1], m[2] .+ s[2] .* [1,-1,-1,1,1], color="black", linestyle=:dashdot,  label="final perturbation")

timesteps = [0.8, 0.5, .25, .1, .05, 0.01]
for (i,dt) = enumerate(timesteps)
    v = deviation(zeros(2); dt)
    n = length(timesteps)
    color = RGB(i/n, 1-i/n, 1-i/n)
    plot!(v[1,:], v[2,:], label="Euler, Δt = $dt"; color)
end

title!("Lokta Volterra Trajectories")
imgdir = "$(@__DIR__)/imgs"
savefig("$imgdir/trajectory.png")
display(Plots.current())

#=
#if q₀
plsot(u_ref.u[1,:], u_ref.u[2,:])

err(m) = abs.(m.-q_ref)./abs(q_ref)

solver = Euler()
timesteps = [1e-1, 1e-2, 1e-3, 1e-4]
for dt in timesteps

end
=#
#=
_curves = [[], []]; labels = []

msizes = [1e-1, 1e-2, 1e-3 ]
nsamples = [50, 100, 200, 400, 800]
for dt in msizes
    q = []; t = []
    for n = nsamples
        ys = zeros(10); τs = zeros(10)
        for i=1:10
            x = sample_perturbation(n); x = vcat(x, -x)
            τ = @elapsed y = mean(quantity(x, dt))
            ys[i] = y; τs[i] = τ
        end
        append!(q, [ys])
        append!(t, [mean(τs)])
    end
    append!(labels, ["Antithetic MC, h=$(dt)"])
    append!(_curves[1], [q])
    append!(_curves[2], [t])
end

plot(); xlims!(extrema(nsamples))
for (ys,l) in zip(_curves[1], labels)
    errorline!(nsamples, err.(ys), label="$l")
end
xlabel!("number of samples")
yaxis!("error")
xaxis!(:log)
display(Plots.current())


plot(); xlims!(extrema(nsamples))
for (t,l) in zip(_curves[2], labels)
    plot!(nsamples, t, label="$l")
end
xlabel!("number of samples")
xlabel!("computation time")
xaxis!(:log)
display(Plots.current())
=#
