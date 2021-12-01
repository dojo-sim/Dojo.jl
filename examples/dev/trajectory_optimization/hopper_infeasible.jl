# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
using IterativeLQR

# Open visualizer
vis = Visualizer()
open(vis)

include(joinpath(module_dir(), "examples", "loader.jl"))

# System
gravity = -9.81
Δt = 0.05
mech = gethopper(Δt = Δt, g = gravity, contact = true, damper = 1.0)
initializehopper!(mech)

## state space
n = minCoordDim(mech)
m = 3 + n

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

z1 = max2min(mech, hopper_offset_state(0.0, 0.0, 0.0))
zM = max2min(mech, hopper_offset_state(0.5, 0.5, 0.5))
zT = max2min(mech, hopper_offset_state(0.5, 0.5, 0.0))

u_control = [0.0; 0.0; mech.g * mech.Δt; zeros(n)]
u_mask = [0 0 0 1 0 0 0;
		  0 0 0 0 1 0 0;
		  0 0 0 0 0 0 1]

# Set random seed
Random.seed!(0)

# Model
function fd(y, x, u, w)
    u_control = u[1:3]
    s = u[3 .+ (1:n)]
	z = simon_step!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)
	y .= copy(max2min(mech, z)) + s
end

function fdx(fx, x, u, w)
	u_control = u[1:3]
    s = u[3 .+ (1:n)]
	fx .= copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	u_control = u[1:3]
    s = u[3 .+ (1:n)]
	∇u = copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[2])
	fu .= [∇u * u_mask' I(n)]
end

# Time
T = 21
h = mech.Δt

dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]

# Initial conditions, controls, disturbances
ū = [[0.0; 0.0; mech.g * mech.Δt + 0.0 * randn(1)[1]; zeros(n)] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
storage = generate_storage(mech, [min2max(mech, x) for x in x̄])
visualize(mech, storage; vis = vis)

# Objective
ot1 = (x, u, w) -> transpose(x - zM) *
	Diagonal(Δt * [0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) *
	(x - zM) + transpose(u) * Diagonal(Δt * [0.01; 0.01; 0.01; 1*ones(n)]) * u
ot2 = (x, u, w) -> transpose(x - zT) *
	Diagonal(Δt * [0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) *
	(x - zT) + transpose(u) * Diagonal(Δt * [0.01; 0.01; 0.01; 1*ones(n)]) * u
oT = (x, u, w) -> transpose(x - zT) *
	Diagonal(Δt * [0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) *
	(x - zT)

ct1 = Cost(ot1, n, m, d)
ct2 = Cost(ot2, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]

# Constraints
function goal(x, u, w)
    Δ = x - zT
    return [Δ[collect(1:6)]; Δ[collect(12 .+ (1:2))]]
end
# Constraints Slack
function slack(x, u, w)
    Δ = u[3 .+ (1:n)]
    return Δ
end

cont = Constraint(slack, n, m)
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)


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
    ρ_scale=3.0)

x_sol, u_sol = IterativeLQR.get_trajectory(prob)
storage = generate_storage(mech, [min2max(mech, x) for x in x_sol])
visualize(mech, storage, vis = vis)

norm([norm(u[3 .+ (1:n)], Inf) for u in u_sol], Inf)

storage_ref = deepcopy(storage)
fxx = prob.m_data.model_deriv.fx[1]

rank(fxx)
indi = [1,2,3,4,5,6,11,12,13]
indi = [1,2,3,8,9,10,11,12,13]
ind = [indi; 13 .+ indi]
cond(fxx)
cond(fxx[indi, indi])
cond(fxx[ind, ind])

plot(Gray.(fxx[indi, indi]))
plot(Gray.(fxx[14:15, 14:15]))
plot(Gray.(fxx))

fxx[14:15, 14:15]
