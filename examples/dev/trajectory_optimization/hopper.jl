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
Δt = 0.05

# Parameters
leg_axis = [0.0; 0.0; 1.0]
leg_length_nominal = 0.5
body_radius = 0.1
foot_radius = 0.05
body_mass = 1.0 
foot_mass = 0.1 

# Links
origin = Origin{Float64}()
body = Sphere(body_radius, body_mass)
foot = Sphere(foot_radius, foot_mass)
links = [body, foot]

# Joint Constraints
joint_origin_body = EqualityConstraint(Floating(origin, body))
joint_body_foot = EqualityConstraint(Prismatic(body, foot, leg_axis; p1=szeros(Float64, 3), p2=szeros(Float64, 3), damper=0.1) )
eqcs = [joint_origin_body, joint_body_foot]

# Contact 
contact_normal = [0.0; 0.0; 1.0]
friction_coefficient = 0.5
contineqcs = contactconstraint(foot, contact_normal, friction_coefficient, p = [0.0; 0.0; 0.0])

# Mechanism
mech = Mechanism(origin, links, eqcs, [contineqcs], g=gravity, Δt=Δt)

# origin to body
setPosition!(mech.origin, mech.bodies[3], Δx=[0.0; 0.0; leg_length_nominal])
setVelocity!(mech.bodies[3], v = [0.0; 0.0; 0.0], ω = [0.0; 0.0; 0.0])
mech.bodies[3].state.xc
mech.bodies[3].state.vc
mech.bodies[3].state.qc
mech.bodies[3].state.ωc

# body to foot
setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; -leg_length_nominal], Δq=UnitQuaternion(RotX(0.0)))
setVelocity!(mech.bodies[4], v = zeros(3), ω = zeros(3))
mech.bodies[4].state.xc
mech.bodies[4].state.vc
mech.bodies[4].state.qc
mech.bodies[4].state.ωc

# controller 
function controller!(mech, k)
    j1 = geteqconstraint(mech, mech.eqconstraints[1].id)
    j2 = geteqconstraint(mech, mech.eqconstraints[2].id)

    setForce!(mech, j1, [0.0; 0.0; 0.0; 0.0; 0.0; 0.0])
    # setForce!(mech, j2, SA[(k < 5 ? 2.5 : (k < 8 ? -1.5 : 0.0)) * mech.g * mech.Δt])
    setForce!(mech, j2, SA[mech.g * mech.Δt])

    return
end 

# simulate
storage = simulate!(mech, 20 * mech.Δt, controller!, record = true, verbose=true, solver = :mehrotra!)#, btol=grad_btol, undercut=grad_undercut)

# visualize
visualize(mech, storage, vis = vis)

# plot(hcat(storage.q[1]...)', linetype=:steppost)

## state space 
n = 13 * 2 
m = 3

function step1!(mech::Mechanism, z, u; btol=1.0e-6, undercut=Inf)
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
        j2 = geteqconstraint(mech, mech.eqconstraints[2].id)
    
        setForce!(mech, j1, [0.0; 0.0; 0.0; u[1]; u[2]; 0.0])
        setForce!(mech, j2, SA[u[3]])
    
        return
    end 

    # simulate  
    storage = simulate!(mech, mech.Δt, 
        controller!, record=true, verbose=false, solver=:mehrotra!, btol=btol, undercut=undercut)
    
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

    return nextstate
end

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

# mid state 
x_shift = 0.5
y_shift = 0.5
z_shift = 0.5

xMb1 = x1b1 + [x_shift; y_shift; z_shift]
vMb1 = [0.0; 0.0; 0.0] 
qMb1 = [1.0; 0.0; 0.0; 0.0]
ωMb1 = [0.0; 0.0; 0.0]
zMb1 = [xMb1; vMb1; qMb1; ωMb1]

xMb2 = x1b2 + [x_shift; 0.0; z_shift]
vMb2 = [0.0; 0.0; 0.0] 
qMb2 = [1.0; 0.0; 0.0; 0.0]
ωMb2 = [0.0; 0.0; 0.0]
zMb2 = [xMb2; vMb2; qMb2; ωMb2]

zM = [zMb1; zMb2] 

# target state 
x_shift = 0.5
y_shift = 0.5
z_shift = 0.0

xTb1 = x1b1 + [x_shift; y_shift; z_shift]
vTb1 = [0.0; 0.0; 0.0] 
qTb1 = [1.0; 0.0; 0.0; 0.0]
ωTb1 = [0.0; 0.0; 0.0]
zTb1 = [xTb1; vTb1; qTb1; ωTb1]

xTb2 = x1b2 + [x_shift; 0.0; z_shift]
vTb2 = [0.0; 0.0; 0.0] 
qTb2 = [1.0; 0.0; 0.0; 0.0]
ωTb2 = [0.0; 0.0; 0.0]
zTb2 = [xTb2; vTb2; qTb2; ωTb2]

zT = [zTb1; zTb2]

u_control = [0.0; 0.0; mech.g * mech.Δt]
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
struct HopperMax{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int
    mech
end

eval_btol = 1.0e-4
eval_undercut = Inf

function fd(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	return step1!(model.mech, x, u, btol=eval_btol, undercut=eval_undercut)
end

grad_btol = 1.0e-3
grad_undercut = 1.5

function fdx(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:(end-model.m)], w[(end-model.m+1):end], btol=grad_btol, undercut=grad_undercut), [x; u])[:, 1:(end-model.m)]
end

function fdu(model::HopperMax{Midpoint, FixedTime}, x, u, w, h, t)
	return fdjac(w -> step1!(mech, w[1:(end-model.m)], w[(end-model.m+1):end], btol=grad_btol, undercut=grad_undercut), [x; u])[:, (end-model.m+1):end]
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
	analytical_dynamics_derivatives = true);

# Solve
@time stats = constrained_ddp_solve!(prob,
    verbose = true,
    grad_tol = 1.0e-3,
	max_iter = 100,
    max_al_iter = 2,
	ρ_init = 1.0,
    ρ_scale = 10.0,
	con_tol = 1.0e-3)

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





