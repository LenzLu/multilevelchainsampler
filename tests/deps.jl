using Pkg
using Test
pkgs = """
    IterTools
    LinearAlgebra    
    Statistics    
    Graphs
    OrdinaryDiffEq
    NonlinearSolve
    QuasiMonteCarlo
    Makie 
    WGLMakie
    GraphMakie
"""
pkgs = filter(x -> x != "" && ! occursin("#", x), 
    split(replace(pkgs, " " => ""),"\n"))

# Activate
project = abspath(dirname(@__DIR__)) 
Pkg.activate(project)

# Install all packages
for pkg in pkgs Pkg.add(pkg) end

# Precompile / test if packages are installed successfully 
@testset "Packages" begin 
    for pkg in pkgs
        @test ( eval(
            quote using $(Symbol(pkg)) end
        ); true) 
    end
end;