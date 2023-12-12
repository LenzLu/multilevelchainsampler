function analyse_nodal_stability(
    g::PowerGrid, node::Int64=1, n=256; t_final=120.0, tolerance=0.1,
    interactive = true)
    N = nv(g.grid)

    δx = sample_perturbations(1, n)
    println("Calculating basin stability ")
    @time S = basin_stability(g, [node], δx; t_final, tolerance)
    
    fig = Figure(aspect=0.3); 

    ax_graph = Axis(fig[1,1])
    graphplot!(ax_graph, g.grid, nlabels=map(string, 1:N))

    ax_pert = Axis(fig[1,2])
         
    U =.! S 
    δϕ, δω = δx
    scatter!(ax_pert, δϕ[S], δω[S], color=:green, label="stable")
    scatter!(ax_pert, δϕ[U], δω[U], color=:tomato, label="unstable")
    axislegend(ax_pert)
        
    s = mean(S)
    ax_pert.xlabel = L"Phase perturbation $δϕ$"
    ax_pert.ylabel = L"Velocity perturbation $δω$"
    ax_pert.title = "Stability of node $node / $N "
    ax_pert.subtitle = L"$ S_1(%$node) \approx %$s \quad (n = %$n) $ "
    

    if interactive
        deregister_interaction!(ax_pert, :rectanglezoom)
        
        ax_traj = Axis(fig[1,3])
        
        
        function plot_trajectories(perturbation)
            empty!(ax_traj); 
            
            ax_traj.xlabel = L"Time $t$"
            ax_traj.ylabel = L"Frequency $ω(t)"
            ax_traj.subtitle = L"$δϕ = %$(round(perturbation[node],digits=5))$ , $δω = %$(round(perturbation[node+1],digits=5)) $" 
            
            xlims!(ax_traj, [-.01 * t_final, t_final] )
            #ylims!(ax_traj, [-10 * tolerance, 10*tolerance])
            lines!(ax_traj,  [0, t_final], repeat([tolerance],2), linestyle=:dash, color=:black)
            lines!(ax_traj,  [0, t_final], repeat([-tolerance],2), linestyle=:dash, color=:black)
                    
            u0 = g.syncstate .+ perturbation
            problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
            solution = solve(problem, Rodas5(); adaptive=true, abstol=1e-6, reltol=1e-6)
            
            #solution = solve(problem, RK4(); adaptive=true, abstol=0.1, reltol=0.1)
            x = [solution.t ... ]
            for i = 1:N
                y = [solution.u[t][i+N] for t=1:length(x) ]
                l = lines!(ax_traj, x, y, label=L"$\omega_{%$i}$")
            end
        end

        plot_trajectories(zeros(2N))
    
        marker = scatter!(ax_pert, [Point2f(0,0)], marker=:xcross, color=:black)
        register_interaction!(ax_pert, :click) do event::MouseEvent, axis   
            if event.type == MouseEventTypes.leftclick 
                pos = mouseposition(axis)
                marker[1][] = [pos]

                perturbation = zeros(2N)
                perturbation[[node, node+N]]   .= pos
                plot_trajectories(perturbation)
                
                sleep(0.01)
            end
        end
    
    end
    return fig
end

function analyse_network(g, interactive=true)
    fig = Figure(aspect=0.6); 
    ax = Axis(fig[1,1])
    
    empty!(ax)
    graphplot!(ax, g.grid, nlabels=map(string, 1:N))
    
    if interactive 
        ax_traj = Axis(fig[1,2])
        
        empty!(ax_traj)
        ax_traj.xlabel = L"Time $t$"
        ax_traj.ylabel = L"Frequency $ω(t)"
        
        t_final = 120.0; tolerance = 0.05
        xlims!(ax_traj, [-.01 * t_final, t_final] )
        ylims!(ax_traj, [-10 * tolerance, 10*tolerance])
        lines!(ax_traj,  [0, t_final], repeat([tolerance],2), linestyle=:dot, color=:black)
        lines!(ax_traj,  [0, t_final], repeat([-tolerance],2), linestyle=:dot, color=:black)
                
        u0 = g.syncstate .+ rand(2N) 
        problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
        solution = solve(problem, Rodas5(); adaptive=true, abstol=1e-6, reltol=1e-6)
        
        #solution = solve(problem, RK4(); adaptive=true, abstol=0.1, reltol=0.1)
        x = [solution.t ... ]
        for i = 1:N
            y = [solution.u[t][i+N] for t=1:length(x) ]
            l = lines!(ax_traj, x, y, label=L"$\omega_{%$i}$")
        end
    end
    return fig
end