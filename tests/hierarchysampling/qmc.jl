using CairoMakie 
using QuasiMonteCarlo

function qmc_points(N, sampler)
    x = QuasiMonteCarlo.sample(N, 2, sampler)
    p = [ Point2f(x[1,i], x[2,i]) for i=1:N ]
    return p
end

function plot_samples(sampler, name)

    outdir = "$(@__DIR__)/outputs"
    points = Observable(Point2f[])
    fig, ax, col = scatter(points, axis=(; limits=(0,1,0,1) ))
    record(fig, "$outdir/$name.mp4", 2:1025) do frame
        p = qmc_points(frame, sampler)
        points[] = p
        N = length(p)
        level = convert(Int, floor(log2(N)))
        notify(points)
        ax.title = "N = $N ≤ 2^$level"
    end
end

plot_samples(SobolSample(), "sobol")
plot_samples(HaltonSample(), "halton")
