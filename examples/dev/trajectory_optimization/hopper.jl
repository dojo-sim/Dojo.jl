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
include(joinpath(module_dir(), "examples", "dev", "trajectory_optimization", "utils.jl"))

# System
gravity = -9.81
Δt = 0.05
mech = gethopper(Δt=Δt, g=gravity)
initializehopper!(mech)

## state space
n = 13 * length(mech.bodies)
m = 3

function hopper_initial_state()
    # initial state
    x1b1 = [0.0; 0.0; 0.5]
    v1b1 = [0.0; 0.0; 0.0]
    q1b1 = [1.0; 0.0; 0.0; 0.0]
    ω1b1 = [0.0; 0.0; 0.0]
    z1b1 = [x1b1; v1b1; q1b1; ω1b1]

    x1b2 = [0.0; 0.0; 0.0]
    v1b2 = [0.0; 0.0; 0.0]
    q1b2 = [1.0; 0.0; 0.0; 0.0]
    ω1b2 = [0.0; 0.0; 0.0]
    z1b2 = [x1b2; v1b2; q1b2; ω1b2]

    z1 = [z1b1; z1b2]
end

function hopper_offset_state(x_shift, y_shift, z_shift)
    z = hopper_initial_state()
    shift = [x_shift; y_shift; z_shift]
    z[1:3] += shift
    z[13 .+ (1:3)] += shift

    return z
end

z1 = hopper_initial_state()
zM = hopper_offset_state(0.5, 0.5, 0.5)
zT = hopper_offset_state(0.5, 0.5, 0.0)

u_control = [0.0; 0.0; mech.g * mech.Δt]
z = [copy(z1)]
for t = 1:5
    znext = step!(mech, z[end], u_control, control_inputs=hopper_inputs!)
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
struct HopperMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech
end

# eval_btol = 1.0e-4
# eval_undercut = Inf

function fd(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	return step!(model.mech, x, u, control_inputs=hopper_inputs!)
end

# grad_btol = 1.0e-3
# grad_undercut = 1.5

function fdx(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	# return fdjac(w -> step!(mech, w[1:(end-model.m)], w[(end-model.m+1):end],
    #     btol=grad_btol, undercut=grad_undercut, control_inputs=hopper_inputs!),
    #     [x; u])[:, 1:(end-model.m)]
    step_grad_x!(model.mech, x, u, control_inputs=hopper_inputs!)
end

function fdu(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	# return fdjac(w -> step!(mech, w[1:(end-model.m)], w[(end-model.m+1):end],
    #     btol=grad_btol, undercut=grad_undercut, control_inputs=hopper_inputs!),
    #     [x; u])[:, (end-model.m+1):end]
    return step_grad_u!(model.mech, x, u, control_inputs=hopper_inputs!)
end

n, m, d = 26, 3, 0
model = HopperMax{Midpoint, FixedTime}(n, m, d, mech);

# Time
T = 21
h = mech.Δt

# Initial conditions, controls, disturbances
ū = [[0.0; 0.0; mech.g * mech.Δt + 0.0 * randn(1)[1]] for t = 1:T-1]
w = [zeros(model.d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w, h, T)
step!(mech, z1, ū[1], control_inputs=hopper_inputs!)
step_grad_x!(mech, z1, ū[1], control_inputs=hopper_inputs!)
step_grad_u!(mech, z1, ū[1], control_inputs=hopper_inputs!)


# Objective
qt1 = [0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
qt2 = [0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
Q = [(t < T ? h * Diagonal([qt1; qt2])
        : h * Diagonal([qt1; qt2])) for t = 1:T]
q = [-2.0 * Q[t] * (t < 11 ? zM : zT) for t = 1:T]

R = [h * Diagonal([0.1; 0.1; 0.01]) for t = 1:T-1]
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

# Constraints
p = [t < T ? 0 : 12 for t = 1:T]
info_t = Dict()
info_T = Dict(:xT => zT)
con_set = [StageConstraint(p[t], t < T ? info_t : info_T) for t = 1:T]

function c!(c, cons::StageConstraints, x, u, t)
	T = cons.T
	p = cons.con[t].p

	if t == T
        Δ = x - zT
		c[1:6] = Δ[1:6]
        c[6 .+ (1:6)] = Δ[13 .+ (1:6)]
	end
end

prob = problem_data(model, obj, con_set, copy(x̄), copy(ū), w, h, T,
    n=[n for t = 1:T], m = [m for t = 1:T-1],
	analytical_dynamics_derivatives = true);
prob.m_data;
prob.m_data.dyn_deriv.fu[4]

# Solve
@time stats = constrained_ddp_solve!(prob,
    verbose = true,
    grad_tol = 1.0e-3,
	max_iter = 100,
    max_al_iter = 2,
	ρ_init = 1.0,
    ρ_scale = 10.0,
	con_tol = 1.0e-3)

x̄, ū = nominal_trajectory(prob)

storage = generate_storage(mech, x̄)
visualize(mech, storage, vis = vis)
