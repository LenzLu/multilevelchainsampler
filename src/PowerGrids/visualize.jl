function analyse_nodal_stability(
    g::PowerGrid, node=1, n=256; t_final=120.0, tolerance=0.1,
    interactive = true)

    δx = sample_perturbations(1, n)
    S = basin_stability(g, [1], δx; t_final, tolerance)
    U = .! S

    fig = Figure(aspect=0.6); ax = Axis(fig[1,1])
    
    δϕ, δω = δx; N = nv(g.grid)
    scatter!(ax, δϕ[S], δω[S], color=:green, label="stable")
    scatter!(ax, δϕ[U], δω[U], color=:tomato, label="unstable")
    axislegend(ax)
    
    s = mean(S)
    ax.xlabel = L"Phase perturbation $δϕ$"
    ax.ylabel = L"Velocity perturbation $δω$"
    ax.title = "Stability of node $node"
    ax.subtitle = L"$ S_1 \approx %$s \quad (n = %$n) $ "

    if interactive
        deregister_interaction!(ax, :rectanglezoom)
        
        function plot_perturbation(perturbation)
            empty!(ax_traj)
            a = g.damping[node]; K = round(mean( g.coupling ), digits=3)
            #text!( ax_traj, L"$\alpha = %$a, K = %$K $", pos=(0.0, 0.9*tolerance), space=:data )
            ax_traj.xlabel = L"Time $t$"
            ax_traj.ylabel = L"Frequency $ω(t)$"
            ax_traj.subtitle = L"$δϕ = %$(round(perturbation[1],digits=5))$ , $δω = %$(round(perturbation[2],digits=5)) $" 


            xlims!(ax_traj, [-.01 * t_final, t_final] )
            ylims!(ax_traj, [-10 * tolerance, 10*tolerance])
            lines!(ax_traj,  [0, t_final], repeat([tolerance],2), linestyle=:dash, color=:black)
            lines!(ax_traj,  [0, t_final], repeat([-tolerance],2), linestyle=:dash, color=:black)
                    
            u0 = g.syncstate .+ perturbation
            problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
            solution = solve(problem, RK4(); adaptive=true, abstol=0.1, reltol=0.1)
            x = [solution.t ... ]
            for i = 1:N
                y = [solution.u[t][i+N] for t=1:length(x) ]
                l = lines!(ax_traj, x, y, label=L"$\omega_{%$i}$")
            end
        end

        ax_traj = Axis(fig[1,2])
        plot_perturbation(zeros(2N))
    
        marker = scatter!(ax, [Point2f(0,0)], marker=:xcross, color=:black)
        register_interaction!(ax, :click) do event::MouseEvent, axis   
            if event.type == MouseEventTypes.leftclick 
                pos = mouseposition(axis)
                marker[1][] = [pos]

                perturbation = zeros(2N)
                perturbation[[node, node+N]]   .= pos
                plot_perturbation(perturbation)
                
                sleep(0.01)
            end
        end
    
    end
    return fig
end