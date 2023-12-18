function plot_state!(ax::Axis, g::PowerGrid; ax_colorbar=nothing)
    node_attr = (; color=g.power, colormap=:bluesreds )
    W = abs.(g.power); Wmin, Wmax = extrema(W)
    node_size = 10 .+ ( abs.(W .- Wmin)) ./ (Wmax .- Wmin) .* 20
    plt = graphplot!(ax, g.grid, nlabels=map(string, 1:nv(g.grid)); node_attr, node_size)
    return plt
end 