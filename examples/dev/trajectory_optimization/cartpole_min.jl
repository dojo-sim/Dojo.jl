# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end
# module_dir
# # Activate package
# using Pkg
# Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

using IterativeLQR

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "dev", "trajectory_optimization", "utils.jl"))

# System
gravity = -9.81
Δt = 0.1
mech = getcartpole(Δt=Δt, g=gravity)
initializecartpole!(mech)

## state space
# n = 13 * 2
n = 4
m = 1
d = 0

function cartpole_initial_state(;pendulum_length=1.0)
    # initial state
    x1b1 = [0.0; 0.0; 0.0]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; -0.5 * pendulum_length]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function cartpole_goal_state(;pendulum_length=1.0)
    # target state
    xTb1 = [0.0; 0.0; 0.0]
    vTb1 = [0.0; 0.0; 0.0]
    qTb1 = [1.0; 0.0; 0.0; 0.0]
    ωTb1 = [0.0; 0.0; 0.0]
    zTb1 = [xTb1; vTb1; qTb1; ωTb1]

    xTb2 = [0.0; 0.0; 0.5 * pendulum_length]
    vTb2 = [0.0; 0.0; 0.0]
    qTb2 = [0.0; 1.0; 0.0; 0.0]
    ωTb2 = [0.0; 0.0; 0.0]
    zTb2 = [xTb2; vTb2; qTb2; ωTb2]

    zT = [zTb1; zTb2]
end

z1 = max2min(mech, cartpole_initial_state())
zT = max2min(mech, cartpole_goal_state())

u_control = [0.0]
u_mask = [1 0]

function fd(y, x, u, w)
	z = simon_step!(mech, min2max(mech,x), u_mask'*u, ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)
	y .= copy(max2min(mech, z))
end

function fdx(fx, x, u, w)
	fx .= copy(getMinGradients!(mech, min2max(mech,x), u_mask'*u, ϵ = 1e-5, btol = 1e-3, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	∇u = copy(getMinGradients!(mech, min2max(mech,x), u_mask'*u, ϵ = 1e-5, btol = 1e-3, undercut = 1.5, verbose = false)[2])
	fu .= ∇u * u_mask'
end

# Time
T = 26
h = mech.Δt

dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]


# Initial conditions, controls, disturbances
ū = [t < 5 ? 1.0 * rand(m) : (t < 10 ? -1.0 * rand(m) : zeros(m)) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w)
storage = generate_storage(mech, [min2max(mech,x) for x in x̄])
visualize(mech, storage; vis = vis)

# Objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(h * ones(n)) * (x - zT) + transpose(u) * Diagonal(Δt * [0.1]) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(100.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

prob = problem_data(model, obj)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# Solve
IterativeLQR.ilqr_solve!(prob,
	verbose = true,
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=500)

x_sol, u_sol = get_trajectory(prob)
storage = generate_storage(mech, [min2max(mech,x) for x in x_sol])
visualize(mech, storage, vis = vis)

fxx = prob.m_data.model_deriv.fx[1]

rank(fxx)
indi = [1,2,3,4,5,6,11,12,13]
indi = [1,2,3]
ind = [indi; 13 .+ indi]
cond(fxx)
cond(fxx[indi, indi])
cond(fxx[ind, ind])

fxx[indi, indi]
