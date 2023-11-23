using Graphs
using Statistics
using Plots
using GraphRecipes

using Revise
using MultilevelChainSampler
imgdir = "$(@__DIR__)/imgs"

Ks = [0.0:0.1:4.0 ... ]; a = 1.0
C2 = complete_graph(2); W_C2 = [-1,1]
C10 = complete_graph(10); W_C10 = [-ones(5)..., ones(5)...]

S = zeros(2,length(Ks))
@time for (i,K) = enumerate(Ks)
    g = PowerGrid(C2, W_C2, K, a)
    S[1,i] = basin_stability(g; N=512)

    g = PowerGrid(C10, W_C10, K, a)
    S[2,i] = basin_stability(g; N=512)
end
plot()
plot!(Ks, S[1,:], label="\$C_2 \$")
plot!(Ks, S[2,:], label="\$C_{10} \$")

xlabel!("coupling \$K\$")
ylabel!("stability \$S\$")
title!("Stability vs. Coupling, α=$a")
savefig("$imgdir/parameters_K_vs_S.png")
display(Plots.current())
