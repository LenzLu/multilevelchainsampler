include("abstract.jl")
include("dynamics.jl")
include("stability.jl")

using Plots
using Graphs

function plot_solution(sol::ODESolution)
    T = sol.t
    d = length(sol[1])
    N = d÷2

    colors = [RGB(0, i/N, 1 - i/N) for i=1:N ]

    plot(); xlims!(extrema(T))
    title!("Phases")
    for i=1:N
        x = [sol.u[t][i] for t=1:length(T)]
        plot!(T, x, label="phase $i", c=colors[i])
    end
    display(Plots.current())

    plot(); xlims!(extrema(T))
    title!("Velocities")
    for i=1:N
        x = [sol.u[t][i+N] for t=1:length(T)]
        plot!(T, x, label="velocity $i", c=colors[i])
    end
    return Plots.current()

end

function draw_oscillators(G::Graph, ϕ::Vector, W::Vector; layout=spring_layout)
    x,y = layout(G)
    draw_oscillators(G, ϕ, W, x, y)
end

function draw_oscillators(G::Graph, ϕ::Vector, W::Vector, x::Vector, y::Vector)
    # Create plot
    plot()
    xlims!(extrema(x).*1.2)
    ylims!(extrema(y).*1.2)

    # Draw edges
    for e in edges(G)
        u,v = e.src, e.dst
        plot!([x[u], x[v]], [y[u], y[v]], label="", color="gray", alpha=.2)
    end

    # Draw oscillators
    R = min(- -(extrema(x)...), - -(extrema(y)...) )/20.
    t = collect(0:.05:1) .* 2π
    z = hcat(cos.(t) , sin.(t))
    for i=1:nv(G)
        p = [ x[i], y[i] ]
        C = reshape(p, (1, 2)) .+ z .* R     # circle
        color =  if (W[i] > 0) "green" else "blue" end
        plot!(C[:,1], C[:,2], label="", color=color, alpha=.4, lw=2)
        if W[i] > 0
            C = reshape(p, (1, 2)) .+ z .* .9*R    # circle
            plot!(C[:,1], C[:,2], label="", color=color, alpha=.4)
        end
        u = [ cos(ϕ[i]), sin(ϕ[i]) ]
        v = p .+ u .* R  # phase
        plot!([p[1],v[1]], [p[2],v[2]], label="", color="red")
        annotate!(p[1], p[2], "$i", 10, color=color)
    end
    return Plots.current()
end

function plot_steady(grid::EnergyGrid)
    u_fix = steady_state(grid)
    ϕ = u_fix[1:nv(grid.G)]
    draw_oscillators(grid.G, ϕ, grid.W)
end
