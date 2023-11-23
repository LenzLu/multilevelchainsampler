using Graphs
using OrdinaryDiffEq
using Statistics
using Plots

using Revise
using MultilevelChainSampler
imgdir = "$(@__DIR__)/imgs"

function mse(f::Function, y_true; n_trials=100)
   square_errors = [
        ( f() .- y_true ).^2
        for k=1:n_trials
   ]
   return mean(square_errors)
end

N=24; K = 1.0; a = 0.5
G = complete_graph(N)
W = ones(nv(G)); W[1:nv(G)÷2] .= -1
g = PowerGrid(G,W,K,a)
timespan = (0.0, 20.0)

u_fix = synchronous_state(g)

perturbe  = MultilevelChainSampler.sample_perturbation
deviation = MultilevelChainSampler.nodal_deviation

println("Calculating reference stability")
perturbation = perturbe(10000)
dev = deviation(g,1; u_fix, perturbation,
      solver=Rodas5(),
      reltol=1e-4, abstol=1e-5)
timed = @elapsed S = mean( dev .< 0.05)
println("Finished after $timed seconds.")
println("S_ref ", S, "\n")

nsamples = [1:10:100 ... ]
errors = zeros(4, length(nsamples))

dev1 = (n) -> deviation(g, 1,perturbe(n); u_fix,
   solver=RK4(), abstol=0.75, reltol=0.75)
dev2 = (n) -> deviation(g, 1,perturbe(n); u_fix,
      solver=RK4(), abstol=0.1, reltol=0.1)
dev3 = (n) -> deviation(g, 1,perturbe(n); u_fix,
   solver=RK4(), abstol=0.05, reltol=0.05)
@time for (i,n) = enumerate(nsamples) # calculate MSE
   println("n = $n")
   stab = () -> mean( dev1(n) .< 0.05 )
   errors[1, i] = mse(stab, S)
   stab = () -> mean( dev2(n) .< 0.05 )
   errors[2, i] = mse(stab, S)
   stab = () -> mean( dev3(n) .< 0.05 )
   errors[3, i] = mse(stab, S)
end

## Plot
plot()
plot!(nsamples, errors[1,:], label="RK4, tol=0.75")
plot!(nsamples, errors[2,:], label="RK4, tol=0.1")
plot!(nsamples, errors[3,:], label="RK4, tol=0.05")
plot!(nsamples, 0.2 .* 1.0 ./ nsamples, label="\$O(N^{-1})\$", color="black", linestyle=:dash )
xlabel!("Number of samples \$N\$")
ylabel!("Mean squared error (MSE)")
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
title!("Estimator error")
savefig("$imgdir/stability_convergence.png")
display(Plots.current())
