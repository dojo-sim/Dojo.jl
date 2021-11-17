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
include(joinpath(module_dir(), "examples", "loader.jl"))

# System 
gravity = -9.81 
Δt = 0.1

# Parameters
slider_axis = [0.0; 1.0; 0.0]
pendulum_axis = [1.0; 0.0; 0.0]
slider_length = 1.0
pendulum_length = 1.0
width, depth, height = 0.1, 0.1, 0.1
slider_mass = 1.0 
pendulum_mass = 1.0 

# Links
origin = Origin{Float64}()
slider = Box(width, slider_length, height, slider_mass)
pendulum = Box(width, depth, pendulum_length, pendulum_mass)
links = [slider, pendulum]

# Joint Constraints
joint_origin_slider = EqualityConstraint(Prismatic(origin, slider, slider_axis; p1=szeros(Float64, 3), p2=szeros(Float64, 3)))
joint_slider_pendulum = EqualityConstraint(Revolute(slider, pendulum, pendulum_axis; p1=szeros(Float64, 3), p2=[0.0; 0.0; 0.5 * pendulum_length]))
eqcs = [joint_origin_slider, joint_slider_pendulum]

# Mechanism
mech = Mechanism(origin, links, eqcs, g=gravity, Δt=Δt)

# origin to slider
setPosition!(mech.origin, mech.bodies[3])
setVelocity!(mech.bodies[3], v = [0.0; 0.0; 0.0], ω = zeros(3))
mech.bodies[3].state.xc
mech.bodies[3].state.vc
mech.bodies[3].state.qc
mech.bodies[3].state.ωc

# slider to pendulum
setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; -0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
setVelocity!(mech.bodies[4], v = zeros(3), ω = zeros(3))
mech.bodies[4].state.xc
mech.bodies[4].state.vc
mech.bodies[4].state.qc
mech.bodies[4].state.ωc

# setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; 0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
# setVelocity!(mech.bodies[4], v = zeros(3), ω = zeros(3))

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
m = 1#isempty(mech.eqconstraints) ? 0 : sum(getcontroldim.(mech.eqconstraints))

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

       # controller 
    function controller!(mech, k)
        j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
        # j2 = geteqconstraint(mech, mech.eqconstraints[2].id)

        setForce!(mech, j1, SVector{1}(u[1]))
        # setForce!(mech, j2, SA[u2])

        return
    end 

    # simulate  
    storage = simulate!(mech, mech.Δt, 
        controller!, record=true, verbose=false, solver=:mehrotra!)
    
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

u_control = 0.0
z = [copy(z1)]
for t = 1:5
    znext = step1!(mech, z[end], [u_control]) 
    push!(z, znext)
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
struct CartpoleMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech
end

model = CartpoleMax{Midpoint, FixedTime}(1,1,1,mech);

function fd(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
	return step1!(model.mech, x, u)
end

function fdx(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(model.mech, w[1:end-1], w[end:end]), [x; u])[:, 1:end-1]
end

function fdu(model::CartpoleMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(model.mech, w[1:end-1], w[end:end]), [x; u])[:, end:end]
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

# Solve
@time ddp_solve!(prob,
    max_iter = 100, verbose = true, linesearch = :armijo)

x, u = current_trajectory(prob)
x̄, ū = nominal_trajectory(prob)

function generate_storage(x) 
    steps = length(x) 
    nbodies = 2
    storage = Storage{Float64}(steps, nbodies)

    for t = 1:steps 
        off = 0
        for (i, body) in enumerate(mech.bodies)
            storage.x[i][t] = x[t][off .+ (1:3)]
            storage.v[i][t] = x[t][off .+ (4:6)]
            storage.q[i][t] = UnitQuaternion(x[t][off .+ (7:10)]...)
            storage.ω[i][t] = x[t][off .+ (11:13)]
            off += 13
        end
    end

    return storage
end

storage = generate_storage(x̄) 
@show x̄[end]

visualize(mech, storage, vis = vis)







