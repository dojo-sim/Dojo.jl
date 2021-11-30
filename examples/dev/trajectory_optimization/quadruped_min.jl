# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "dev", "trajectory_optimization", "utils.jl"))

using IterativeLQR

# System
gravity = -9.81
Δt = 0.05
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, damper = 10.0, spring = 0.0)
initialize!(mech, :quadruped, tran = [0.00,0,0.], v = [0.5,0,0.])

ineqcs = collect(mech.ineqconstraints)

contact_location(mech, ineqcs[1])
contact_location(mech, ineqcs[2])
contact_location(mech, ineqcs[3])
contact_location(mech, ineqcs[4])

@elapsed storage = simulate!(mech, 0.05, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)



################################################################################
# Quadruped Trajectory Generator
################################################################################

# Trajectory template
function trunk_trajectory(x::T, z::T, r::T, N::Int) where T
	# x = forward displacement
	# z = vertical displacement
	# r = radius of the foot traj
	# trunk trajectory decomposed into N steps (N+1 waypoints)
	X = Vector(range(-r + x, r + x, length = N))
	trunk_traj = [[x, 0.0, z] for x in X]
	return trunk_traj
end

function low_foot_trajectory(r::T, x::T, y::T, N::Int) where T
	# r = radius of the foot traj
	# x = forward displacement
	# y = lateral displacement
	# low foot trajectory decomposed into N steps (N+1 waypoints)
	X = Vector(range(-r + x, r + x, length = N+1))
	low_traj = [[x, y, 0.0] for x in X]
	return low_traj
end

function high_foot_trajectory(r::T, x::T, y::T, N::Int) where T
	# high foot trajectory decomposed into N steps (N+1 waypoints)
	α = range(0, π, length = N+1)
	θ = π * (1 .- cos.(α)) / 2
	X = [x + r*cos(θi) for θi in θ]
	Z = [r*sin(θi)/2 for θi in θ] # the facto 2 flatten the foo trajectory
	high_traj = [[X[i], y, Z[i]] for i = 1:N+1]
	return high_traj
end

function foot_trajectory(r::T, x::T, y::T, N::Int) where T
	low_traj = low_foot_trajectory(r, x, y, N)
	high_traj = high_foot_trajectory(r, x, y, N)
	traj = [low_traj[1:end-1]; high_traj[1:end-1]]
	return traj
end

# Inverse kinematics: given the position of the trunk, finds the knee and hip angles that will put the foot
# at the p_foot location
function inverse_kinematics(mechanism::Mechanism, p_trunk, p_foot; leg::Symbol = :FR)
	# starting point of the local search
	θ = [0.95, -1.5*0.95] # θhip, θknee
	for k = 1:10
		err = IKerror(mechanism, p_trunk, p_foot, θ; leg = leg)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(θ -> IKerror(mechanism, p_trunk, p_foot, θ; leg = leg), θ)
		θ -= ∇ \ err
	end
	return θ
end

function IKerror(mechanism::Mechanism, p_trunk, p_foot, θ; leg::Symbol = :FR)
	setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [p_trunk; zeros(3)])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_thigh_joint"), [θ[1]])
	setPosition!(mechanism, geteqconstraint(mechanism, String(leg)*"_calf_joint"), [θ[2]])

	foot = getbody(mechanism, String(leg)*"_calf")
	ineqcs = collect(mechanism.ineqconstraints)
	ineqc = ineqcs[findfirst(x -> x.parentid == foot.id, ineqcs)]
	p = contact_location(ineqc, foot)
	err = p - p_foot
	return err[[1,3]]
end

function quadruped_trajectory(mechanism::Mechanism{T}, r, x, z; N = 6, Ncycles = 1) where T
	pFR = [+ 0.13, - 0.13205, 0.]
	pFL = [+ 0.13, + 0.13205, 0.]
	pRR = [- 0.23, - 0.13205, 0.]
	pRL = [- 0.23, + 0.13205, 0.]

	t = reverse(foot_trajectory(r, 0.0, 0.0, N))
	t_dephased = [t[N+1:2N]; t[1:N]]
	# Foot positions
	tFR = [pFR + ti for ti in t]
	tFL = [pFL + ti for ti in t_dephased]
	tRR = [pRR + ti for ti in t_dephased]
	tRL = [pRL + ti for ti in t]

	# Leg angles
	p_trunk = [0,0,z]
	θFR = [inverse_kinematics(mechanism, p_trunk, tFR[i], leg = :FR) for i = 1:2N]
	θFL = [inverse_kinematics(mechanism, p_trunk, tFL[i], leg = :FL) for i = 1:2N]
	θRR = [inverse_kinematics(mechanism, p_trunk, tRR[i], leg = :RR) for i = 1:2N]
	θRL = [inverse_kinematics(mechanism, p_trunk, tRL[i], leg = :RL) for i = 1:2N]

	# Minimal Coordinates
	X = [Vector{T}([p_trunk; zeros(3); zeros(6);
		zeros(2); θFR[i][1]; 0; θFR[i][2]; 0;
		zeros(2); θFL[i][1]; 0; θFL[i][2]; 0;
		zeros(2); θRR[i][1]; 0; θRR[i][2]; 0;
		zeros(2); θRL[i][1]; 0; θRL[i][2]; 0;
		]) for i = 1:2N]

	X = vcat([deepcopy(X) for i = 1:Ncycles]...)
	# adding forward velocity
	for i = 1:2N*Ncycles
		X[i][1] += (i-1) * 2r / N
	end
	return X
end

################################################################################
# Compute trajectory
################################################################################

# discretization
N = 6 # half period of the 1-step loop

# quadruped geometry
y = 0.1 # lateral displacement

# trunk trajectory
z = 0.29 # trunk height
v = 0.40 # forward speed

# Half disk foot trajectory
r = 0.10 # foot traj radius


X = quadruped_trajectory(mech, r, 0.0, z; N = N, Ncycles = 1)
storage = generate_storage(mech, [min2max(mech, x) for x in X])
visualize(mech, storage, vis = vis)


collect(mech.ineqconstraints)
p_trunk = [0,0,0.31]
p_foot = [0.2,0,0]
inverse_kinematics(mech, p_trunk, p_foot)



N = 6
low_traj = low_foot_trajectory(r, 0.0, y, Nfoot)
high_traj = high_foot_trajectory(r, 0.0, y, Nfoot)
traj = foot_trajectory(r, 0.0, y, Nfoot)
trunk_traj = trunk_trajectory(0.0, z, r, Nfoot)

plot()
scatter!([p[1] for p in low_traj], [p[3] for p in low_traj], )
scatter!([p[1] for p in high_traj], [p[3] for p in high_traj], )
scatter!([p[1] for p in traj], [p[3] for p in traj], )



plot()
scatter!([p[1] for p in tFR], [p[3] for p in tFR])
scatter!([p[1] for p in tFL], [p[3] for p in tFL])
scatter!([p[1] for p in tRR], [p[3] for p in tRR])
scatter!([p[1] for p in tRR], [p[3] for p in tRR])





## state space
Nb = length(mech.bodies)
n = minCoordDim(mech)
m = 12

function potato_dynamics(x, u, Δt, m, g)
	# x = [x,y,z,ẋ,ẏ,ż]
	gv = [0, 0, g]
	ẋ = [x[4:6]; u ./ (m*Δt) + gv]
	x̄ = x + ẋ * Δt
	return x̄
end

trunk = getbody(mech, "trunk")
x2_trunk = trunk.state.x2[1]
v15_trunk = trunk.state.v15


U_potato = [fill([0,0,8], 6); fill(zeros(3), 15)]
x_potato = [x2_trunk; v15_trunk]
X_potato = []
for t = 1:21
	mass = sum(getfield.(mech.bodies, :m))
	alt = x_potato[3] - 0.48
	if t > 13
		u_potato = -[0, 0, mech.Δt * mass * mech.g + 20*alt + 10*x_potato[6]]
	else
		u_potato = U_potato[t]
	end
	push!(X_potato, x_potato)
	x_potato = potato_dynamics(x_potato, u_potato, mech.Δt, mass, mech.g)
end
plot()
plot!([x[1] for x in X_potato], linewidth = 5.0)
plot!([x[2] for x in X_potato], linewidth = 5.0)
plot!([x[3] for x in X_potato], linewidth = 5.0)


zref = []
for t = 1:21
	initialize!(mech, :quadruped, tran = X_potato[t][1:3] - [0,0,0.31], v = X_potato[t][4:6], θ = 0.95)
	push!(zref, getMaxState(mech))
end
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)
zref = [max2min(mech, z) for z in zref]

initialize!(mech, :quadruped, tran = [0.00,0,0.0], v = [0.5,0,0], θ = 0.95)
z1 = max2min(mech, getMaxState(mech))
visualizeMaxCoord(mech, min2max(mech, z1), vis)

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

mech = getmechanism(:quadruped, Δt = Δt, g = gravity, damper = 1000.0, spring = 30.0)
initialize!(mech, :quadruped)
@elapsed storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
ugc = gravity_compensation(mech)

mech = getmechanism(:quadruped, Δt = Δt, g = gravity, damper = 10.0, spring = 0.0)

u_control = ugc[6 .+ (1:12)]
u_mask = [zeros(12,6) I(m)]

z = [copy(z1)]
for t = 1:5
    znext = max2min(mech, simon_step!(mech, min2max(mech, z[end]), u_mask'*u_control))
    push!(z, znext)
end

# Model
function fd(y, x, u, w)
	z = simon_step!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)
	y .= copy(max2min(mech, z))
end

function fdx(fx, x, u, w)
	fx .= copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	∇u = copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[2])
	fu .= ∇u * u_mask'
end


# Time
T = 21
h = mech.Δt

n, m, d = minCoordDim(mech), 12, 0
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]


# Initial conditions, controls, disturbances
ū = [u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Rollout
initialize!(mech, :quadruped, tran = [0.00,0,0.0], v = [0.0,0,0], θ = 1.1)
z1 = max2min(mech, getMaxState(mech))

x̄ = rollout(model, z1, ū, w)
# simon_step!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-6, undercut = 1.5, verbose = false)
# getGradients!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-3, undercut = 1.5, verbose = false)
storage = generate_storage(mech, [min2max(mech, x) for x in x̄])
visualize(mech, storage; vis = vis)

# Objective
# qt1 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
# qt2 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
body_scale = [1; 0.1ones(12)]
qt = vcat([body_scale[i] * [0.1 * ones(3); 0.001 * ones(3); 0.1 * ones(4); 0.01 * ones(3)] for i = 1:Nb]...)
qt = [0.1 * ones(3); 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(24)]

# ot1 = (x, u, w) -> transpose(x - zM) * Diagonal(Δt * qt) * (x - zM) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
# ot2 = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
# oT = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT)
ots = [(x, u, w) -> transpose(x - zref[t]) * Diagonal(Δt * qt) * (x - zref[t]) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u for t = 1:20]
oT = (x, u, w) -> transpose(x - zref[end]) * Diagonal(Δt * qt) * (x - zref[end])

# ct1 = Cost(ot1, n, m, d)
# ct2 = Cost(ot2, n, m, d)
# cT = Cost(oT, n, 0, 0)
cts = Cost.(ots, n, m, d)
cT = Cost(oT, n, 0, 0)
# obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]
obj = [cts..., cT]

# Constraints
function goal(x, u, w)
	# Δ = x - zT
    Δ = x - zref[end]
    return Δ[collect(1:3)]
end

cont = Constraint()
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# Solve
IterativeLQR.constrained_ilqr_solve!(prob,
    verbose = true,
	linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0)

x_sol, u_sol = get_trajectory(prob)
storage = generate_storage(mech, [min2max(mech, x) for x in x_sol])
visualize(mech, storage, vis = vis)
