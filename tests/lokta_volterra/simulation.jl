

dt = 1e-3
δu = sample_perturbation(100)
n = size(δu,1)
V = [deviation(δu[i,:]; dt) for i=1:n]
Q = [v[1,end] for v in V]
