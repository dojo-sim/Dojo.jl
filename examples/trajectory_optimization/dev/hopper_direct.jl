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

using DirectTrajectoryOptimization

include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "src", "optional_components", "trajopt_utils.jl"))

# System
gravity = -9.81
Δt = 0.05
mech = getraiberthopper(Δt = Δt, g = gravity, damper=0.0)
initializeraiberthopper!(mech)

## state space
n = minCoordDim(mech)
m = 3

function raiberthopper_initial_state()
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

function raiberthopper_offset_state(x_shift, y_shift, z_shift)
    z = raiberthopper_initial_state()
    shift = [x_shift; y_shift; z_shift]
    z[1:3] += shift
    z[13 .+ (1:3)] += shift
    return z
end

z1 = max2min(mech, raiberthopper_offset_state(0.0, 0.0, 0.0))
zM = max2min(mech, raiberthopper_offset_state(0.25, 0.0, 0.5))
zT = max2min(mech, raiberthopper_offset_state(0.5, 0.0, 0.0))

# z1 = max2min(mech, raiberthopper_offset_state(0.0, 0.0, 0.0))
# zM = max2min(mech, raiberthopper_offset_state(0.0, 0.0, 0.0))
# zT = max2min(mech, raiberthopper_offset_state(0.0, 0.0, 0.0))

u_control = [0.0; 0.0; mech.g * mech.Δt]
u_mask = [0 0 0 1 0 0 0;
		  0 0 0 0 1 0 0;
		  0 0 0 0 0 0 1]

# Set random seed
Random.seed!(0)

# Model
nx, nu, nw = minCoordDim(mech), 3, 0
nu += nx
function f(d, y, x, u, w)
    u_ctrl = u[1:3]
    s = u[3 .+ (1:nx)]
    z = step!(mech, min2max(mech, x), u_mask'*u_ctrl, ϵ = 1e-6, btol = 3e-4, undercut = 1.5, verbose = false)
	d .= y - max2min(mech, z) + s
end

function fz(dz, y, x, u, w)
    u_ctrl = u[1:3]
    s = u[3 .+ (1:nx)]
	dx, du = getMinGradients!(mech, min2max(mech, x), u_mask'*u_ctrl, ϵ = 1e-6, btol = 3e-3, undercut = 1.5, verbose = false)
    dz .= [-dx -du * transpose(u_mask) I(nx) I(nx)]
end

# Time
T = 21
h = mech.Δt


# ## model
dt = Dynamics(f, fz, nx, nx, nu)
dyn = [dt for t = 1:T-1]
model = DynamicsModel(dyn)

# ## objective
ot1 = (x, u, w) -> 10.0 * transpose(x - zM) * Diagonal([1.0; 1.0; 100.0; 10.0; 10.0; 10.0; 1.0e-1; 1.0e-1; 1.0e-1; 10.0; 100.0; 10.0; 100.0; 0.1]) * (x - zM) + 1.0 * transpose(u[1:3]) * Diagonal([10.0; 10.0; 1.0]) * u[1:3] + 1.0 * transpose(u[3 .+ (1:nx)]) * u[3 .+ (1:nx)]
ot2 = (x, u, w) -> 1.0e-1 * transpose(x - zT) * Diagonal([1.0; 1.0; 1.0; 10.0; 10.0; 10.0; 1.0e-1; 1.0e-1; 1.0e-1; 10.0; 100.0; 10.0; 1.0; 1.0]) * (x - zT) + 1.0 * transpose(u[1:3]) * Diagonal([10.0; 10.0; 1.0]) * u[1:3] + 1.0 * transpose(u[3 .+ (1:nx)]) * u[3 .+ (1:nx)]
oT = (x, u, w) -> 10.0 * transpose(x - zT) * Diagonal([10.0; 10.0; 10.0; 10.0; 10.0; 10.0; 1.0e-1; 1.0e-1; 1.0e-1; 10.0; 100.0; 10.0; 10.0; 1.0]) * (x - zT)

# ot1 = (x, u, w) -> transpose(x - zM) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zM) + transpose(u[1:3]) * 100.0 * Diagonal([0.01; 0.01; 0.01]) * u[1:3] + 1.0 * transpose(u[3 .+ (1:nx)]) * u[3 .+ (1:nx)]
# ot2 = (x, u, w) -> transpose(x - zT) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zT) + transpose(u[1:3]) * 100.0 * Diagonal([0.01; 0.01; 0.01]) * u[1:3] + 1.0 * transpose(u[3 .+ (1:nx)]) * u[3 .+ (1:nx)]
# oT = (x, u, w) -> transpose(x - zT) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zT)

ct1 = Cost(ot1, nx, nu, nw, [t for t = 1:10])
ct2 = Cost(ot2, nx, nu, nw, [t for t = 11:20])
cT = Cost(oT, nx, 0, 0, [T])
obj = [ct1, ct2, cT]

# ## constraints
function _x_init(x, u, w)
    x - z1
end
function slack(x, u, w)
    u[3 .+ (1:nx)]
end
function _x_goal(x, u, w)
    (x - zT)[1:3]
end
con1 = StageConstraint(_x_init, nx, nu, nw, [1], :equality)
cont = StageConstraint(slack, nx, nu, nw, [t for t = 1:T-1], :equality)
conT = StageConstraint(_x_goal, nx, 0, 0, [T], :equality)
cons = ConstraintSet([con1, cont])#, conT])

# ## problem
trajopt = TrajectoryOptimizationProblem(obj, model, cons)
s = Solver(trajopt, options=Options(
    tol=1.0e-2,
    constr_viol_tol=1.0e-2,
    compl_inf_tol=1.0e-2,
    acceptable_tol=1.0e-2,
    acceptable_constr_viol_tol = 1.0e-2,
    acceptable_compl_inf_tol = 1.0e-2,
    mu_target = 1.0e-2
))

# ## initialize
x_interpolation = [linear_interpolation(z1, zM, 11)..., linear_interpolation(zM, zT, 11)[2:end]...]
ū = [[0.0; 0.0; mech.g * mech.Δt + 0.0 * randn(1)[1]; 1.0e-3 * randn(nx)] for t = 1:T-1]
z0 = zeros(s.p.num_var)
for (t, idx) in enumerate(s.p.trajopt.model.idx.x)
    z0[idx] = x_interpolation[t]
end
for (t, idx) in enumerate(s.p.trajopt.model.idx.u)
    z0[idx] = ū[t]
end
DirectTrajectoryOptimization.initialize!(s, z0)

# ## solve
@time solve!(s)

# ## solution
@show trajopt.x[1]
@show trajopt.x[T]

@show sum([sum(abs.(u[3 .+ (1:nx)])) for u in trajopt.u[1:end-1]])
# ## state
# plot(hcat(trajopt.x...)')

# # ## control
# plot(hcat(trajopt.u[1:end-1]..., trajopt.u[end-1])', linetype = :steppost)

storage = generate_storage(mech, [min2max(mech, x) for x in trajopt.x])
visualize(mech, storage, vis = vis)
