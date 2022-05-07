using Test
using LinearAlgebra
using Symbolics
using BenchmarkTools

include("layer.jl")
include("net.jl")

################################################################################
# Layer
################################################################################
ni = 3
no = 6
nθ = no * (ni + 1)
activation0(x) = tanh.(x)

xi0 = rand(ni)
θ0 = rand(nθ)

layer0 = Layer148(ni, no, activation0)
layer0.xi = xi0
layer0.θ = θ0

xo0 = evaluation(ni, no, xi0, θ0, activation0)
evaluation!(layer0)
@benchmark $evaluation!($layer0)
@test norm(xo0 - layer0.xo, Inf) < 1e-5

J0 = FiniteDiff.finite_difference_jacobian(xi0 -> evaluation(ni, no, xi0, θ0, activation0), xi0)
jacobian_input!(layer0)
@benchmark $jacobian_input!($layer0)
@test norm(J0 - layer0.∂xo∂xi, Inf) < 1e-5

J0 = FiniteDiff.finite_difference_jacobian(θ0 -> evaluation(ni, no, xi0, θ0, activation0), θ0)
jacobian_parameters!(layer0)
@benchmark $jacobian_parameters!($layer0)
@test norm(J0 - layer0.∂xo∂θ, Inf) < 1e-5


################################################################################
# Net
################################################################################
function local_evaluation(net, xi, θ)
    evaluation!(net, xi, θ)
    return copy(get_output(net))
end

ni = 3
n0 = 15
xi0 = rand(ni)
net0 = Net148(ni, no, dim_layers=[6,9,12], activations=[x->x, x->x, x->x, x->x])
θ0 = rand(parameter_dimension(net0))

jacobian_input!(net0, xi0, θ0)
J0 = FiniteDiff.finite_difference_jacobian(xi0 -> local_evaluation(net0, xi0, θ0), xi0)
@test norm(net0.∂xo∂xi - J0, Inf) < 1e-5


jacobian_parameters!(net0, xi0, θ0)
J0 = FiniteDiff.finite_difference_jacobian(θ0 -> local_evaluation(net0, xi0, θ0), θ0)
@test norm(net0.∂xo∂θ - J0, Inf) < 1e-5
