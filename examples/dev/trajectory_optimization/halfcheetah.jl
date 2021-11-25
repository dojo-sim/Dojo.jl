# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end
# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat
using Colors
using GeometryBasics
using Rotations
using Parameters
using Symbolics
using Random
using LinearAlgebra
using ForwardDiff

# Open visualizer
vis = Visualizer()
open(vis)

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

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "dev", "trajectory_optimization", "utils.jl"))


# System
gravity = -9.81
Δt = 0.05
mech = getmechanism(:halfcheetah, Δt = Δt, g = gravity)
initialize!(mech, :halfcheetah)

## state space
Nb = length(mech.bodies)
n = 13 * Nb
m = 3

z1 = halfcheetahState(x = 0.0, z = 0.00, θ = 0.0)
zM = halfcheetahState(x = 0.0, z = 0.40, θ = 0.0)
zT = halfcheetahState(x = 0.5, z = 0.00, θ = 0.0)


function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = controldim(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.eqconstraints
        nu = controldim(eqc)
        if eqc.parentid != nothing
            body = getbody(mechanism, eqc.parentid)
            rot = eqc.constraints[2]
            A = Matrix(nullspacemat(rot))
            Fτ = springforce(mechanism, eqc, body)
            F = Fτ[1:3]
            τ = Fτ[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end

mech = getmechanism(:halfcheetah, Δt = Δt, g = gravity, damper = 100.0, spring = 1000.0)
initialize!(mech, :halfcheetah, x = 0.0, z = 0.0, θ = 0.0)
@elapsed storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
ugc = gravity_compensation(mech)

mech = getmechanism(:halfcheetah, Δt = Δt, g = gravity, damper = 10.0, spring = 1000.0)


u_control = ugc
u_mask = I(length(u_control))

z = [copy(z1)]
for t = 1:5
    znext = simon_step!(mech, z[end], u_mask'*u_control)
    push!(z, znext)
end


# Set random seed
Random.seed!(0)
# Model
struct HopperMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech
end

function fd(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	return simon_step!(model.mech, x, u_mask'*u, ϵ = 1e-5, btol = 1e-5, undercut = 1.5, verbose = false)
end

function fdx(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	∇x, ∇u = getGradients!(model.mech, x, u_mask'*u, ϵ = 1e-4, btol = 1e-3, undercut = 1.5, verbose = false)
	return ∇x
end

function fdu(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	∇x, ∇u = getGradients!(model.mech, x, u_mask'*u, ϵ = 1e-4, btol = 1e-3, undercut = 1.5, verbose = false)
	return ∇u * u_mask'
end

n, m, d = 13Nb, 9, 0
model = HopperMax{Midpoint, FixedTime}(n, m, d, mech);

# Time
T = 21
h = mech.Δt

# Initial conditions, controls, disturbances
ū = [u_control for t = 1:T-1]
w = [zeros(model.d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w, h, T)
simon_step!(model.mech, z1, u_mask'*u_control, ϵ = 1e-6, btol = 1e-6, undercut = 1.5, verbose = false)
getGradients!(model.mech, z1, u_mask'*u_control, ϵ = 1e-6, btol = 1e-3, undercut = 1.5, verbose = false)
storage = generate_storage(mech, x̄)
visualize(mech, storage; vis = vis)

# Objective
qt1 = [0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
qt2 = [0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
body_scale = [1; 0.2ones(6)]
qt = vcat([body_scale[i] * [0.1 * ones(3); 0.001 * ones(3); 0.1 * ones(4); 0.01 * ones(3)] for i = 1:Nb]...)
Q = [(t < T ? h * Diagonal(qt)
        : h * Diagonal(qt)) for t = 1:T]
q = [-2.0 * Q[t] * (t < 11 ? zM : zT) for t = 1:T]

# R = [h * Diagonal([0.1; 0.1; 0.01]) for t = 1:T-1]
R = [h * Diagonal(0.01 * ones(length(u_control))) for t = 1:T-1]
r = [-2.0 * R[t] * u_control for t = 1:T-1]
# r = [-0.0 * R[t] * u_control for t = 1:T-1]

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
stats = constrained_ddp_solve!(prob,
    verbose = true,
    grad_tol = 1.0e-3,
	max_iter = 100,
    max_al_iter = 5,
	ρ_init = 1.0,
    ρ_scale = 10.0,
	con_tol = 1.0e-3)

x̄, ū = nominal_trajectory(prob)

storage = generate_storage(mech, x̄)
visualize(mech, storage, vis = vis)


minimalCoordinates(mechanism, eqc)
