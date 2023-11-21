using Test

@testset "trajectory"  begin
    include("deviation/trajectory.jl")
end

@testset "pushforward" begin
    include("deviation/pushforward.jl")
end

@testset "rates" begin
    include("convergence/rates.jl")
end
