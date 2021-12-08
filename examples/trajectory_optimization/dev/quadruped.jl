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
include(joinpath(module_dir(), "src", "optional_components", "trajopt_utils.jl"))

using IterativeLQR

# System
gravity = -9.81
Δt = 0.05
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, damper = 1.0, spring = 100.0)
initialize!(mech, :quadruped, v = [0.5, 0, 0])
# @elapsed storage = simulate!(mech, 0.5, record = true, solver = :mehrotra!, verbose = false)
# visualize(mech, storage, vis = vis)

## state space
Nb = length(mech.bodies)
n = 13 * Nb
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
	alt = x_potato[3] - 0.40
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
	initialize!(mech, :quadruped, tran = X_potato[t][1:3] - [0,0,0.23], v = X_potato[t][4:6], θ = 0.95)
	push!(zref, getMaxState(mech))
end
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)


initialize!(mech, :quadruped, tran = [0.00,0,0.0], v = [0.5,0,0], θ = 0.95)
z1 = getMaxState(mech)
visualizeMaxCoord(mech, z1, vis)

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

mech = getmechanism(:quadruped, Δt = Δt, g = gravity, damper = 10.0, spring = 30.0)

u_control = ugc[6 .+ (1:12)]
u_mask = [zeros(12,6) I(m)]

z = [copy(z1)]
for t = 1:5
    znext = step!(mech, z[end], u_mask'*u_control)
    push!(z, znext)
end

# Set random seed
Random.seed!(0)

# Model
function fd(y, x, u, w)
	y .= copy(step!(mech, x, u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false))
end

function fdx(fx, x, u, w)
	fx .= copy(getMaxGradients!(mech, x, u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	∇u = copy(getMaxGradients!(mech, x, u_mask'*u, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[2])
	fu .= ∇u * u_mask'
end

# Time
T = 21
h = mech.Δt

n, m, d = 13Nb, 12, 0
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]


# Initial conditions, controls, disturbances
ū = [u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Rollout
x̄ = rollout(model, z1, ū, w)
# step!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-6, undercut = 1.5, verbose = false)
# getGradients!(model.mech, x, u_mask'*u_control, ϵ = 1e-6, btol = 1e-3, undercut = 1.5, verbose = false)
storage = generate_storage(mech, x̄)
visualize(mech, storage; vis = vis)

# Objective
# qt1 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
# qt2 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
body_scale = [1; 0.1ones(12)]
qt = vcat([body_scale[i] * [0.1 * ones(3); 0.001 * ones(3); 0.1 * ones(4); 0.01 * ones(3)] for i = 1:Nb]...)

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
	Δ = x - zT
    Δ = x - zref[end]
    return Δ[collect(1:6)]
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
storage = generate_storage(mech, x_sol)
visualize(mech, storage, vis = vis)
