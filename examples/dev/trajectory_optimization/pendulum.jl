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

# System 
gravity = -9.81 
Δt = 0.1

# Parameters
pendulum_axis = [1.0; 0.0; 0.0]
pendulum_length = 1.0
width, depth, height = 0.1, 0.1, 0.1
pendulum_mass = 1.0 

# Float type 
T = Float64 

# Links
origin = Origin{T}()
pendulum = Box(width, depth, pendulum_length, pendulum_mass)
links = [pendulum]

# Joint Constraints
joint_slider_pendulum = EqualityConstraint(Revolute(origin, pendulum, pendulum_axis; p1=zeros(3), p2=[0.0; 0.0; 0.5 * pendulum_length]))
eqcs = [joint_slider_pendulum]

# Mechanism
mech = Mechanism(origin, links, eqcs, g=gravity, Δt=Δt)

# origin to pendulum
setPosition!(mech.origin, mech.bodies[2], Δx=[0.0; 0.0; -0.5 * pendulum_length])
setVelocity!(mech.bodies[2], v = [0.0; 0.0; 0.0], ω = zeros(3))

setPosition!(mech.origin, mech.bodies[2], Δx=[0.0; 0.0; 0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
setVelocity!(mech.bodies[2], v = zeros(3), ω = zeros(3))

u_control = 0.0

# controller
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)

    u1 = u_control

    setForce!(mech, j1, [u1])

    return
end 

# simulate
storage = simulate!(mech, 1 * mech.Δt, controller!, record = true, solver = :mehrotra!)

# x2 = mech.bodies[2].state.xc
# v1 = mech.bodies[2].state.vc
# qc = mech.bodies[2].state.qc
# q2 = [qc.w; qc.x; qc.y; qc.z]
# ω1 = mech.bodies[2].state.ωc
# [x2; v1; q2; ω1]

plot(hcat(storage.x[1]...)', color=:black, width=2.0)

# visualize
visualize(mech, storage, vis = vis)

## state space 
n = 13 * length(mech.bodies)
m = isempty(mech.eqconstraints) ? 0 : sum(getcontroldim.(mech.eqconstraints))

# function step!(mech::Mechanism, z::Vector{T}, u::Vector{T}) 
#     # set data
#     data = [z; u] 

#     off = 0
  
#     for body in mech.bodies
#         x2, v15, q2, ω15 = unpackdata(data[off+1:end]); off += 13
#         body.state.xc = x2 - v15 * mech.Δt
#         body.state.vc = v15
#         body.state.qc = UnitQuaternion(q2...) * ωbar(-ω15, mech.Δt) * mech.Δt / 2.0
#         body.state.ωc = ω15
#     end

#     function controller!(mech, k)
#         j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
#         setForce!(mech, j1, SVector{1}(u))
#         return
#     end 

#     # simulate  
#     storage = simulate!(mech, mech.Δt, 
#         controller!, 
#         record=true, solver=:mehrotra!)

#     # next state
#     nextstate = Vector{T}()  
#     for body in mech.bodies
#         x3 = body.state.xk[1]
#         v25 = body.state.vsol[2]
#         _q3 = body.state.qk[1]
#         q3 = [_q3.w; _q3.x; _q3.y; _q3.z]
#         ω25 = body.state.ωsol[2]
#         push!(nextstate, [x3; v25; q3; ω25]...)
#     end

#     return nextstate
# end

function step1!(mech::Mechanism, z, u)
    # set data
    data = [z; u] 

    off = 0
  
    for body in mech.bodies
        x2, v15, q2, ω15 = unpackdata(data[off+1:end]); off += 13
        body.state.xc = x2 - v15 * mech.Δt
        body.state.vc = v15
        body.state.qc = UnitQuaternion(q2...) * ωbar(-ω15, mech.Δt) * mech.Δt / 2.0
        body.state.ωc = ω15
    end

    function controller!(mech, k)
        j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
        setForce!(mech, j1, SVector{1}(u))
        return
    end 

    # simulate  
    storage = simulate!(mech, mech.Δt, 
        controller!, 
        record=true, verbose=false, solver=:mehrotra!)
    
    # next state
    nextstate = []#Vector{Float64}()  

    for body in mech.bodies
        x3 = body.state.xk[1]
        v25 = body.state.vsol[2]
        _q3 = body.state.qk[1]
        q3 = [_q3.w; _q3.x; _q3.y; _q3.z]
        ω25 = body.state.ωsol[2]
        push!(nextstate, [x3; v25; q3; ω25]...)
    end

    # # gradient 
    # data = getdata(deepcopy(mech))
    # sol = getsolution(deepcopy(mech))
    # δ = sensitivities(deepcopy(mech), sol, data)

    return nextstate#, δ
end

# initial state 
x1 = [0.0; 0.0; -0.5 * pendulum_length] 
v1 = [0.0; 0.0; 0.0] 
q1 = [1.0; 0.0; 0.0; 0.0]
ω1 = [0.0; 0.0; 0.0] 
z1 = [x1; v1; q1; ω1]

# target state 
xT = [0.0; 0.0; 0.5 * pendulum_length]
vT = [0.0; 0.0; 0.0] 
qT = [0.0; 1.0; 0.0; 0.0]
ωT = [0.0; 0.0; 0.0]
zT = [xT; vT; qT; ωT]

z = [copy(z1)]
# δz = Matrix{T}[]
for t = 1:5
    # znext, δznext = step!(mech, z[end], [u_control]) 
    znext = step1!(mech, z[end], [u_control]) 
    push!(z, znext)
    # push!(δz, δznext)
end 

# plot!(hcat(z...)[1:3, :]')

# # sensi = -1.0 * (full_matrix(mech.system) \ (full_data_matrix(mech) * transpose(attitudejacobian([z[end-1]; u_control], 1))))[6:11, :]
# data = getdata(deepcopy(mech))
# sol = getsolution(deepcopy(mech))
# δ = sensitivities(deepcopy(mech), sol, data)
# sensi = (δ * transpose(attitudejacobian([z[end-1]; u_control], 1)))[6:11, :]
# # attitudejacobian([z[end-1]; u_control], 1)

# # idx 
# idx_x = collect(1:3) 
# idx_v = collect(4:6) 
# idx_q = collect(7:10) 
# idx_ω = collect(11:13)
# idx_u = collect(14:14) 

# idx_vsol = collect(1:3) 
# idx_ωsol = collect(4:6) 

# # var
# x3 = z[end][idx_x]
# v2 = z[end][idx_v]
# q3 = z[end][idx_q]
# ω2 = z[end][idx_ω]

# # data
# x2 = z[end-1][idx_x]
# v1 = z[end-1][idx_v]
# q2 = z[end-1][idx_q]
# ω1 = z[end-1][idx_ω]

# ∂v2∂x2 = sensi[idx_vsol, idx_x]
# ∂v2∂v1 = sensi[idx_vsol, idx_v]
# ∂v2∂q2 = sensi[idx_vsol, idx_q]
# ∂v2∂ω1 = sensi[idx_vsol, idx_ω]
# ∂v2∂u1 = sensi[idx_vsol, idx_u]

# ∂ω2∂x2 = sensi[idx_ωsol, idx_x]
# ∂ω2∂v1 = sensi[idx_ωsol, idx_v]
# ∂ω2∂q2 = sensi[idx_ωsol, idx_q]
# ∂ω2∂ω1 = sensi[idx_ωsol, idx_ω]
# ∂ω2∂u1 = sensi[idx_ωsol, idx_u]

# ∂x3∂x2 = I(3) + ∂v2∂x2 .* mech.Δt
# ∂x3∂v1 = ∂v2∂v1 .* mech.Δt
# ∂x3∂q2 = ∂v2∂q2 .* mech.Δt
# ∂x3∂ω1 = ∂v2∂ω1 .* mech.Δt
# ∂x3∂u1 = ∂v2∂u1 .* mech.Δt

# ∂q3∂ω2 = Lmat(UnitQuaternion(q2...)) * derivωbar(SVector{3}(ω2), mech.Δt) * mech.Δt / 2

# ∂q3∂x2 = ∂q3∂ω2 * ∂ω2∂x2
# ∂q3∂v1 = ∂q3∂ω2 * ∂ω2∂v1
# ∂q3∂q2 = Rmat(ωbar(ω2, mech.Δt) * mech.Δt / 2) + ∂q3∂ω2 * ∂ω2∂q2
# ∂q3∂ω1 = ∂q3∂ω2 * ∂ω2∂ω1
# ∂q3∂u1 = ∂q3∂ω2 * ∂ω2∂u1

# jac = [∂x3∂x2 ∂x3∂v1 ∂x3∂q2 ∂x3∂ω1 ∂x3∂u1;
#        ∂v2∂x2 ∂v2∂v1 ∂v2∂q2 ∂v2∂ω1 ∂v2∂u1;
#        ∂q3∂x2 ∂q3∂v1 ∂q3∂q2 ∂q3∂ω1 ∂q3∂u1;
#        ∂ω2∂x2 ∂ω2∂v1 ∂ω2∂q2 ∂ω2∂ω1 ∂ω2∂u1]


# norm(z[end] - mystep([z[end-1]; u_control])) < 1.0e-8
# norm((FiniteDiff.finite_difference_jacobian(mystep, [z[end-1]; u_control]) - jac)[7:10, 7:10], Inf)
# norm((FiniteDiff.finite_difference_jacobian(mystep, [z[end-1]; u_control]) - jac)[11:13, 7:10], Inf)
# norm((FiniteDiff.finite_difference_jacobian(mystep, [z[end-1]; u_control]) - jac)[[collect(1:6); collect(11:13)], 7:10], Inf)

# (FiniteDiff.finite_difference_jacobian(mystep, [z[end-1]; u_control]))[7:10, 7:10]
# (jac)[7:10, 7:10]

# (FiniteDiff.finite_difference_jacobian(mystep, [z[1]; u_control]))[1:3, 7:10]
# (jac)[1:3, 7:10]


# Rmat(ωbar(ω2, mech.Δt) * mech.Δt / 2) + ∂q3∂ω2 * ∂ω2∂q2


using Colors
using GeometryBasics
using Rotations
using Parameters
using Symbolics
using Random

## get motion_planning.jl and set path
path_mp = "/home/taylor/Research/motion_planning"
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
struct PendulumMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech::Mechanism
end

function fd(model::PendulumMax{Midpoint, FixedTime}, x, u, w, h, t)
	return step1!(model.mech, x, u)
end

function fdx(model::PendulumMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:end-1], w[end:end]), [x; u])[:, 1:end-1]
end

function fdu(model::PendulumMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:end-1], w[end:end]), [x; u])[:, end:end]
end

n, m, d = 13, 1, 0
model = PendulumMax{Midpoint, FixedTime}(n, m, d, mech)

# Time
T = 11
h = mech.Δt

# Initial conditions, controls, disturbances
ū = [0.0 * rand(model.m) for t = 1:T-1]
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
