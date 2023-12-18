# Random walk Proposal
@kwdef struct CyclicWalker <: AbstractProposalGenerator
    lb = 0.0; ub = 10.0
    step = 1.0
end
function propose(g::CyclicWalker) 
    return g.lb + rand() * (g.ub-g.lb) 
end
function propose(g::CyclicWalker, x)
    y = x + randn() * g.step 
    return g.lb + mod( y - g.lb, g.ub - g.lb )
end
