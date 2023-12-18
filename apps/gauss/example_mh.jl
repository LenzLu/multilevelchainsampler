include("utils.jl")
run(`clear`)
imgdir = "$(@__DIR__)/imgs"

## Sample from MetropolisHastings
println("MetropolisHastings: ")
energy = GaussianEnergy([64])
target = x->evaluate(energy, x) 
simple_energy = EnergyFunction(target, tracked=true)
@time s = sample(simple_energy, MetropolisHastings(100000), CyclicWalker())

println("\tChain length ", length(s))
println("\tRejection rate ", mean(s.rejections))
println("\tMean acceptance ", mean(s[:acceptance]))
println("\tEvaluations ", simple_energy.evaluations)
println("")

fig = Figure()
ax_x, ax_E = plot_chain_statistics!(fig, s)
ax_x.title = "MetropolisHastings"
plot_energy!(ax_x, x->(x .- 5)^2, color=:black, linestyle=:dash, label=L"$E(x)=(x-5)^2$")
plot_energy!(ax_x, target, color=:red, label="target")
Legend(fig[2,1], ax_x, tellwidth = false, tellheight = true, orientation=:horizontal)
save("$imgdir/metropolis_hastings.png", fig)
display(fig)


## Sample from DelayedAcceptance
println("Delayed acceptance: ")
energy = GaussianEnergy([16, 128])
sampler = MetropolisHastings(150000)
@time s = sample(energy, sampler, CyclicWalker())

println("\tChain length ", length(s))
println("\tTotal rejection rate ", mean(s.rejections))
println("\tRejection rate ", mean(hcat(s[:rejections]...), dims=2)[:,1])
println("\tMean acceptance ", mean(s[:acceptance]))
println("\tEvaluations ", energy.evaluations)
println("")

fig = Figure()
ax_x, ax_E = plot_chain_statistics!(fig, s)
ax_x.title = "Delayed Acceptance" 
plot_energy!(ax_x, x->(x .- 5)^2, color=:black, linestyle=:dash, label=L"$E(x)=(x-5)^2$")
target_F = x->evaluate(energy, x, level=2) 
target_C = x->evaluate(energy, x, level=1) 
plot_energy!(ax_x, target_F, color=:red, label="fine")
plot_energy!(ax_x, target_C, color=:cyan, label="course")
Legend(fig[2,1], ax_x, tellwidth = false, tellheight = true, orientation=:horizontal)
save("$imgdir/delayed_acceptance.png", fig)
display(fig)

## Sample from DelayedAcceptance
println("Delayed acceptance (flatter energy):")
target_F = x->evaluate(energy, x, level=2) 
target_C = x->0.3 * evaluate(energy, x, level=1) 
sampler = MetropolisHastings(140000)
e = EnergyFunction(target_C, target_F; tracked=true)
@time s = sample(e, sampler, CyclicWalker())

println("\tChain length ", length(s))
println("\tTotal rejection rate ", mean(s.rejections))
println("\tRejection rate ", mean(hcat(s[:rejections]...), dims=2)[:,1])
println("\tMean acceptance ", mean(s[:acceptance]))
println("\tEvaluations ", e.evaluations)
println("")

fig = Figure()
ax_x, ax_E = plot_chain_statistics!(fig, s)
ax_x.title = "Delayed Acceptance (flatter energy)" 
plot_energy!(ax_x, x->(x .- 5)^2, color=:black, linestyle=:dash, label=L"$E(x)=(x-5)^2$")
plot_energy!(ax_x, target_F, color=:red, label="fine")
plot_energy!(ax_x, target_C, color=:cyan, label="course")
Legend(fig[2,1], ax_x, tellwidth = false, tellheight = true, orientation=:horizontal)
save("$imgdir/delayed_acceptance_flatter.png", fig)
display(fig)

## Sample from DelayedAcceptance
println("Delayed acceptance (even flatter energy):")
energy = GaussianEnergy([32, 64, 128])
target_F = x->evaluate(energy, x, level=3) 
target_M = x->0.3 * evaluate(energy, x, level=2) 
target_C = x->0.1 * evaluate(energy, x, level=1)
sampler = MetropolisHastings(170000)
e = EnergyFunction(target_C, target_M, target_F; tracked=true)
@time s = sample(e, sampler, CyclicWalker())

println("\tChain length ", length(s))
println("\tTotal rejection rate ", mean(s.rejections))
println("\tRejection rate ", mean(hcat(s[:rejections]...), dims=2)[:,1])
println("\tMean acceptance ", mean(s[:acceptance]))
println("\tEvaluations ", e.evaluations)
println("")

fig = Figure()
ax_x, ax_E = plot_chain_statistics!(fig, s)
ax_x.title = "Delayed Acceptance (even flatter energy)" 
plot_energy!(ax_x, x->(x .- 5)^2, color=:black, linestyle=:dash, label=L"$E(x)=(x-5)^2$")
plot_energy!(ax_x, target_F, color=:red, label="fine")
plot_energy!(ax_x, target_M, color=:green, label="medium")
plot_energy!(ax_x, target_C, color=:cyan, label="course")
Legend(fig[2,1], ax_x, tellwidth = false, tellheight = true, orientation=:horizontal)

display(fig)