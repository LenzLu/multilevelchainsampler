include("functionality.jl")

m_ref = quantity(zeros(1,2), 1e-6)[1]
println("True: $m_ref ")

x = sample_perturbation(100); y = quantity(x, 1e-3)
m = mean(y); s = std(y)
println("Vanilla MC: $m, relerr $( abs(m-m_ref)/m_ref ) ")

m = multilevel_estimator(sample_perturbation, quantity)
println("Fixed MLMC: $m, relerr $( abs(m-m_ref)/m_ref )")

m = adaptive_multilevel_estimator(sample_perturbation, quantity; n0=100)
println("Adaptive MLMC: $m, relerr $( abs(m-m_ref)/m_ref )")
