using Dojo
using StaticArrays

# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end
# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(@__DIR__, "..", "..", "loader.jl"))
include(joinpath(module_dir(), "examples/dev/fd_tools.jl"))

# System
gravity = 0.0#-9.81
Δt = 0.1

# Parameters
width, depth, height = 1.0, 1.0, 1.0
mass = 1.0

# Float type
T = Float64

# Links
origin = Origin{T}()
rigidbody = Box(width, depth, height, mass)
links = [rigidbody]

# Joint Constraints
pin_joint = EqualityConstraint(Spherical(origin, rigidbody; p1=zeros(3), p2=[0.0; 0.0; 0.0]))
eqcs = [pin_joint]

# Mechanism
mech = Mechanism(origin, links, eqcs, g=gravity, Δt=Δt)

# origin to pendulum
setPosition!(mech.origin, mech.bodies[2], Δq=UnitQuaternion(RotX(0.2 * π)))
setVelocity!(mech.bodies[2], v = [0.0; 0.0; 0.0], ω = [0.0; 0.0; 0.0])

u_control = [0.0; 0.0; 0.1]

# controller
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)

    setForce!(mech, j1, mech.Δt * u_control)

    return
end

# simulate
storage = simulate!(mech, 10 * mech.Δt, controller!, record = true, solver = :mehrotra!)

# plot(hcat(storage.q[1]...)', color=:black, width=2.0)

# visualize
visualize(mech, storage, vis = vis)

## state space
n = 13 * length(mech.bodies)
m = isempty(mech.eqconstraints) ? 0 : sum(getcontroldim.(mech.eqconstraints))


function step1!(mech::Mechanism, z, u)
    # set data
    data = [z; u]

    off = 0

    for body in mech.bodies
        x2, v15, q2, ω15 = unpackdata(data[off+1:end]); off += 13
        body.state.x1 = x2 - v15 * mech.Δt
        body.state.v15 = v15
        body.state.q1 = UnitQuaternion(q2...) * ωbar(-ω15, mech.Δt) * mech.Δt / 2.0
        body.state.ϕ15 = ω15
    end

    function controller!(mech, k)
        j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
        setForce!(mech, j1, SVector{3}(u))
        return
    end

    # simulate
    storage = simulate!(mech, mech.Δt,
        controller!,
        record=true, verbose=false, solver=:mehrotra!)

    # next state
    nextstate = []#Vector{Float64}()

    for body in mech.bodies
        x3 = body.state.x2[1]
        v25 = body.state.vsol[2]
        _q3 = body.state.q2[1]
        q3 = [_q3.w; _q3.x; _q3.y; _q3.z]
        ω25 = body.state.ϕsol[2]
        push!(nextstate, [x3; v25; q3; ω25]...)
    end

    return nextstate
end

# initial state
x1 = [0.0; 0.0; 0.0]
v1 = [0.0; 0.0; 0.0]
q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 0.0; 0.0]
z1 = [x1; v1; q1; ω1]

# target state
xT = [0.0; 0.0; 0.0]
vT = [0.0; 0.0; 0.0]
_qT = rand(UnitQuaternion)
qT = [_qT.w; _qT.x; _qT.y; _qT.z]
ωT = [0.0; 0.0; 0.0]
zT = [xT; vT; qT; ωT]

z = [copy(z1)]
for t = 1:5
    znext = step1!(mech, z[end], u_control)
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
struct RigidBody{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech::Mechanism
end

function fd(model::RigidBody{Midpoint, FixedTime}, x, u, w, h, t)
	return step1!(model.mech, x, u)
end

function fdx(model::RigidBody{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:(end-model.m)], w[(end-model.m+1):end]), [x; u])[:, 1:(end-model.m)]
end

function fdu(model::RigidBody{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:(end-model.m)], w[(end-model.m+1):end]), [x; u])[:, (end-model.m+1):end]
end

n, m, d = 13, 3, 0
model = RigidBody{Midpoint, FixedTime}(n, m, d, mech)

# Time
T = 11
h = mech.Δt

# Initial conditions, controls, disturbances
ū = [0.01 * rand(model.m) for t = 1:T-1]
w = [zeros(model.d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w, h, T)

# Objective
Q = [(t < T ? Diagonal(1.0e-1 * ones(model.n))
        : Diagonal(1000.0 * ones(model.n))) for t = 1:T]
q = [-2.0 * Q[t] * zT for t = 1:T]

R = [Diagonal(1.0e-4 * ones(model.m)) for t = 1:T-1]
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
    analytical_dynamics_derivatives = true)

# Solve
@time ddp_solve!(prob,
    max_iter = 10, verbose = true, linesearch = :armijo)

x, u = current_trajectory(prob)
x̄, ū = nominal_trajectory(prob)

function generate_storage(x)
    steps = length(x)
    nbodies = 1
    storage = Storage{Float64}(steps, nbodies)

    for t = 1:steps
        storage.x[1][t] = x[t][1:3]
        storage.v[1][t] = x[t][4:6]
        storage.q[1][t] = UnitQuaternion(x[t][7:10]...)
        storage.ω[1][t] = x[t][11:13]
    end

    return storage
end

storage = generate_storage(x)
@show x̄[end]

visualize(mech, storage, vis = vis)
