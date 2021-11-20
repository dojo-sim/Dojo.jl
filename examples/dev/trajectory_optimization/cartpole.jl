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

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "dev", "fd_tools.jl"))
include(joinpath(module_dir(), "examples", "dev", "trajectory_optimization", "utils.jl"))

# System
gravity = -9.81
Δt = 0.1
mech = getcartpole(Δt=Δt, g=gravity)

initializecartpole!(mech)


# controller
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
    j2 = geteqconstraint(mech, mech.eqconstraints[2].id)

    u1 = 0.2
    u2 = 0.0

    setForce!(mech, j1, SA[u1])
    setForce!(mech, j2, SA[u2])

    return
end

# simulate
storage = simulate!(mech, mech.Δt, controller!, record = true, verbose=true, solver = :mehrotra!)

# visualize
visualize(mech, storage, vis = vis)

## state space
n = 13 * 2
m = 1

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

z1 = cartpole_initial_state()
zT = cartpole_goal_state()

u_control = 0.0
z = [copy(z1)]
for t = 1:5
    znext = step!(mech, z[end], [u_control], control_inputs=cartpole_inputs!)
    push!(z, znext)
end

using Colors
using GeometryBasics
using Rotations
using Parameters
using Symbolics
using Random

## get motion_planning.jl and set path
# path_mp = "/home/taylor/Research/motion_planning"
path_mp = joinpath(module_dir(), "..", "motion_planning")
include(joinpath(path_mp, "src/utils.jl"))
include(joinpath(path_mp, "src/time.jl"))
include(joinpath(path_mp, "src/model.jl"))
include(joinpath(path_mp, "src/integration.jl"))
include(joinpath(path_mp, "src/objective.jl"))
include(joinpath(path_mp, "src/constraints.jl"))

# differential dynamic programming
include(joinpath(path_mp, "src/differential_dynamic_programming/ddp.jl"))

using Random, LinearAlgebra, ForwardDiff
Random.seed!(0)

# Model
struct CartpoleMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech
end

function fd(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
	return step!(model.mech, x, u, control_inputs=cartpole_inputs!)
end

function fdx(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
    step_grad_x!(model.mech, x, u, control_inputs=cartpole_inputs!)
end

function fdu(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
    step_grad_u!(model.mech, x, u, control_inputs=cartpole_inputs!)
end

n, m, d = 26, 1, 0

model = CartpoleMax{Midpoint, FixedTime}(n, m, d, mech);

# Time
T = 26
h = mech.Δt

# Initial conditions, controls, disturbances
ū = [t < 5 ? 1.0 * rand(model.m) : (t < 10 ? -1.0 * rand(model.m) : zeros(model.m)) for t = 1:T-1]
w = [zeros(model.d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w, h, T)

# Objective
Q = [(t < T ? Diagonal(h * ones(model.n))
        : Diagonal(1000.0 * ones(model.n))) for t = 1:T]
q = [-2.0 * Q[t] * zT for t = 1:T]

R = [h * Diagonal(1.0e-1 * ones(model.m)) for t = 1:T-1]
r = [zeros(model.m) for t = 1:T-1]

obj = StageCosts([QuadraticCost(Q[t], q[t],
	t < T ? R[t] : nothing, t < T ? r[t] : nothing) for t = 1:T], T)

function g(obj::StageCosts, x, u, t)
	T = obj.T
    if t < T
		Q = obj.cost[t].Q
		q = obj.cost[t].q
	    R = obj.cost[t].R
		r = obj.cost[t].r
        return x' * Q * x + q' * x + u' * R * u + r' * u
    elseif t == T
		Q = obj.cost[T].Q
		q = obj.cost[T].q
        return x' * Q * x + q' * x
    else
        return 0.0
    end
end

# Problem
prob = problem_data(model, obj, copy(x̄), copy(ū), w, h, T,
    analytical_dynamics_derivatives = true);

step!(mech, z1, ū[1], control_inputs=cartpole_inputs!)
step_grad_x!(mech, z1, ū[1], control_inputs=cartpole_inputs!)
step_grad_u!(mech, z1, ū[1], control_inputs=cartpole_inputs!)


# Solve
@time ddp_solve!(prob,
    max_iter = 100, verbose = true, linesearch = :armijo)

x̄, ū = nominal_trajectory(prob)

storage = generate_storage(mech, x̄)

visualize(mech, storage, vis = vis)
