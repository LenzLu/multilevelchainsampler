include("utils.jl")
run(`clear`)

## Sample from WangLandau
println("WangLandau: ")
energy = GaussianEnergy([64])
sampler = WangLandau(0, 100, 10000; nbins=100)
@time s = sample(energy, sampler, CyclicWalker())

println("\tChain length ", length(s))
println("\tRejection rate ", mean(s.rejections))
println("")

fig = Figure()
ax_x, ax_E = plot_chain_statistics!(fig, s; unfolded=false)
ax_x.title = "WangLandau" 
plot_energy!(ax_x, x->(x .- 5)^2, color=:black, linestyle=:dash, label=L"$E = (x-5)^2$")
target = x->evaluate(energy, x) 
plot_energy!(ax_x, target, color=:gray, linestyle=:dash, label=L"E")
Emin,Emax = extrema(s.energies)
lines!(ax_E, [Emin,Emax], repeat([1/(Emax-Emin)],2), color=:blue, linestyle=:dash, label="Flat")
display(fig)



