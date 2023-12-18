function plot_trajectories!(ax::Axis, g::PowerGrid, perturbation::Vector; t_final=120.0, tolerance=0.1)    
    ax.xlabel = L"Time $t$"
    ax.ylabel = L"Frequency $ω(t)"
    
    xlims!(ax, [-.01 * t_final, t_final] )
    #ylims!(ax, [-10 * tolerance, 10*tolerance])
    boundaries = (; linestyle=:dash, color=:black, linewidth=0.5 )
    lines!(ax,  [0, t_final], repeat([tolerance],2); boundaries...)
    lines!(ax,  [0, t_final], repeat([-tolerance],2); boundaries...)
            
    u0 = g.syncstate .+ perturbation
    problem = ODEProblem(swing_dynamics!, u0, (0, t_final), g)
    
    solution = solve(problem, Rodas5(); adaptive=true, abstol=1e-3, reltol=1e-3)
    #solution = solve(problem, RK4(); adaptive=true, abstol=0.2, reltol=0.1)
    
    x = [solution.t ... ]; N = length(u0) ÷ 2
    for i = 1:N
        y = [solution.u[t][i+N] for t=1:length(x) ]
        l = lines!(ax, x, y, label=L"$\omega_{%$i}$")
    end
end
