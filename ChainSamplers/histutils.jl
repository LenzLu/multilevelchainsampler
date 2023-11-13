struct EnergyHistogram
    min::Real
    max::Real
    N::Integer
end

function midpoints(h::EnergyHistogram)
    s = ( h.max - h.min ) / h.N
    m = [h.min+s/2:s:h.max-s/2 ... ]
    return m
end

function which_bin(h::EnergyHistogram, x::Real)
    x = max(h.min, min(h.max, x)) # clamp to bounds
    n = floor((h.N-1) * (x - h.min)/(h.max - h.min))
    n = convert(Integer,n) + 1
    return n
end

function histsample(X,E; n_bins=10)
    H = EnergyHistogram(minimum(E),maximum(E),n_bins)
    bin_idx = [ which_bin(H, e) for e in E]

    indices = -ones(Integer, n_bins)
    for k=1:n_bins
        v = ( bin_idx .== k )
        i = argmax( v )
        if (v[i] > 0)
            indices[k] = i
        end
    end

    return [ (X[i],E[i]) for i in indices if i!=-1 ]
end
