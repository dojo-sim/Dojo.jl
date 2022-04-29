using Pkg
Pkg.activate(joinpath(module_dir(), "examples"))

using Dojo
using IterativeLQR
using RoboDojo
using Plots
using Symbolics
using BenchmarkTools
using LinearAlgebra
using FiniteDiff

const iLQR = IterativeLQR
const RD = RoboDojo

include("methods.jl")


################################################################################
# Dynamics query
################################################################################
dynamics_model = Simulator(RD.halfcheetah, 1, h=h)
dynamics_model.ip.opts.r_tol = 1e-5
dynamics_model.ip.opts.κ_tol = 3e-2
dynamics_model.ip.opts.undercut = 5.0

nq = dynamics_model.model.nq
nx = 2nq
nu = dynamics_model.model.nu
nw = dynamics_model.model.nw



y = zeros(2nq)
dx = zeros(2nq,2nq)
du = zeros(2nq,nu)
x = RD.nominal_state(dynamics_model.model) + 0.1*rand(18)
u = 0.0ones(nu)
w = 0.0ones(nw)

RD.dynamics(dynamics_model, y, x, u, w)
y
RD.dynamics_jacobian_state(dynamics_model, dx, x, u, w)
dx
RD.dynamics_jacobian_input(dynamics_model, du, x, u, w)
du

function explicit_dynamics(dynamics_model, x, u, w)
    y = zeros(2nq)
    RD.dynamics(dynamics_model, y, x, u, w)
    return y
end
dx0 = FiniteDiff.finite_difference_jacobian(x -> explicit_dynamics(dynamics_model, x, u, w), x)
du0 = FiniteDiff.finite_difference_jacobian(u -> explicit_dynamics(dynamics_model, x, u, w), u)
norm(dx - dx0)
norm(du - du0)

# @benchmark $dynamics_jacobian_state($dynamics_model, $dx, $x, $u, $w)
# @benchmark $dynamics_jacobian_control($dynamics_model, $du, $x, $u, $w)
# @benchmark $dynamics($dynamics_model, $y, $x, $u, $w)
