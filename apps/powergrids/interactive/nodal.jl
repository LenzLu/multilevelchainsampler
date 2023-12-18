using WGLMakie
using MultilevelChainSampler

function analyze_nodal_stability(g::PowerGrid, node::Int64=1, n=256; t_final=120.0, tolerance=0.1)
    δx = MultilevelChainSampler.sample_perturbations(1, n)
    println("Calculating stabilities ")

    @time S = basin_stability(g, [node], δx; t_final, tolerance)
    s = mean(S)

    fig = Figure(size=(900, 350))
    
    ax_graph = Axis(fig[1,1])
    plt = plot_state!(ax_graph, g)
    cmap = Makie.extract_colormap(plt.plots[3])
    Colorbar(fig[2,1], colormap=cmap, vertical=false, label="power")

    ax_pert = Axis(fig[1,2])    
    U =.! S 
    δϕ, δω = δx
    scS = scatter!(ax_pert, δϕ[S], δω[S], color=:green, label="stable")
    scU = scatter!(ax_pert, δϕ[U], δω[U], color=:tomato, label="unstable")
    Legend(fig[2,2], ax_pert, tellwidth = false, tellheight = true, orientation=:horizontal)
    ax_pert.xlabel = L"Phase perturbation $δϕ$"
    ax_pert.ylabel = L"Velocity perturbation $δω$"
    ax_pert.title = "Stability of node $node / $N "
    ax_pert.subtitle = L"$ S_1(%$node) \approx %$s \quad (n = %$n) $ "
        
    ax_traj = Axis(fig[1,3])
    plot_trajectories!(ax_traj, g, zeros(2N); t_final, tolerance)
    marker = scatter!(ax_pert, [Point2f(0,0)], marker=:xcross, color=:black)
    ax_traj.title = L"$\delta \phi = 0.0$, $\delta \omega = 0.0 $"
                
    deregister_interaction!(ax_pert, :rectanglezoom)
    register_interaction!(ax_pert, :click) do event::MouseEvent, axis   
        if event.type == MouseEventTypes.leftclick 
            pos = mouseposition(axis)
            marker[1][] = [pos]

            perturbation = zeros(2N)
            perturbation[[node, node+N]] .= pos
            empty!(ax_traj)
            ax_traj.title = L"$δϕ = %$(round(pos[1], digits=3))$, $\delta \omega = %$(round(pos[2], digits=3)) $"
            plot_trajectories!(ax_traj, g, perturbation; t_final, tolerance)
            
            sleep(0.01)
        end
    end 
    return fig
end
