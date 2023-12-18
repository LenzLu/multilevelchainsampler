# sample n perturbations for a group of m nodes
function sample_perturbations(
    m::Int64, n::Int64=128; 
    sampler=HaltonSample(), maxfreq=10.0)
    
    UV = QuasiMonteCarlo.sample(n,2m,SobolSample())
    δϕ = (2 .* UV[1:m    ,:] .- 1 ) .* π
    δω = (2 .* UV[m+1:end,:] .- 1 ) .* maxfreq
    return δx = (δϕ, δω)
end

# apply the perturbations to a group of m nodes
function perturb_nodes(N::Int64, J::Vector{Int64}, δx=sample_perturbation(length(J), 128) )
    m = length(J); @assert all( 1 .<= J .<= N) && J == unique(J)
    δϕ, δω = δx; n = length(δϕ) 
    perturbations = zeros(n, 2N)
    for (i,u)=enumerate(J) 
        perturbations[:, u]   .= δϕ[i, :]
        perturbations[:, u+N] .= δω[i, :]
    end
    return perturbations
end

