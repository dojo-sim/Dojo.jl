using ConstrainedDynamics
using LinearAlgebra

include("examples/chain_in_chain.jl")

ConstrainedDynamics.setentries!(mech)
A1, _, b1 = ConstrainedDynamics.densesystem(mech)
x1 = A1\b1

ConstrainedDynamics.ldu_solve!(mech.system)

A2, x2, b2 = ConstrainedDynamics.densesystem(mech)

@test isapprox(norm(x1-x2), 0.0; atol = 2e-5)