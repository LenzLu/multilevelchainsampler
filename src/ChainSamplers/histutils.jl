struct FixedHistogram{T}
    min::T
    max::T
    N::Int64
end

function midpoints(h::FixedHistogram)
    s = ( h.max - h.min ) / h.N
    m = [h.min+s/2:s:h.max-s/2 ... ]
    return m
end

function which_bin(h::FixedHistogram{T}, x::T) where {T}
    x = max(h.min, min(h.max, x)) # clamp to bounds
    n = floor((h.N-1) * (x - h.min)/(h.max - h.min))
    n = convert(Integer,n) + 1
    return n
end

function histsample(X,E; n_bins=10)
    H = FixedHistogram(minimum(E),maximum(E),n_bins)
    bin_idx = [ which_bin(H, e) for e in E]

    indices = -ones(Int64, n_bins)
    for k=1:n_bins

        v = ( bin_idx .== k )
        i = argmax( v ) # TODO: pick others as well

        if (v[i] > 0)
            indices[k] = i
        end
    end

    valid_samples = (-)
    X,E = [ (X[i],E[i]) for i in indices if (i !=-1) ]
    return X,E
end
