using Test
using Statistics

using MultilevelChainSampler

mutable struct Dummy
    x::Float64 # State
end
struct RandomWalk <: ProposalGenerator{Dummy}
    s::Float64 # Step size
end
MultilevelChainSampler.initialize(R::RandomWalk) = Dummy(randn() * R.s)
MultilevelChainSampler.propose!(R::RandomWalk, d::Dummy) = ( d.x += randn() * R.s; d)

R = RandomWalk(0.1)
chain_length = 100000
sampler = MetropolisHastings(d->0.5 * d.x^2, R, chain_length)
X,E = sample_chain(sampler)
X = [d.x for d in X]
X = X[1000:end]

@test abs( mean(X) ) < 0.1
@test abs( std(X) - 1.0 ) < 0.1
