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
mech = getmechanism(:halfcheetah, Δt = Δt, g = gravity, damper = 1.0, spring = 0.0)
initialize!(mech, :halfcheetah)

## state space
Nb = length(mech.bodies)
n = 13 * Nb
m = 6

z1 = halfcheetahState(x = 0.00, z = 0.00, θ = 0.0)
zM = halfcheetahState(x = 0.25, z = 0.40, θ = 0.0)
zT = halfcheetahState(x = 0.50, z = 0.00, θ = 0.0)


function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = control_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.joints
        nu = control_dimension(eqc)
        if eqc.parentid != nothing
            body = get_body(mechanism, eqc.parentid)
            rot = eqc.constraints[2]
            A = Matrix(nullspace_mask(rot))
            Fτ = apply_spring(mechanism, eqc, body)
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

u_control = ugc[3 .+ (1:6)]
u_mask = [zeros(6,3) I(m)]

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

n, m, d = 13Nb, 6, 0
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
qt1 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
qt2 = [0.1; 0.1; 1.0; 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)]
body_scale = [1; 0.2ones(6)]
qt = vcat([body_scale[i] * [0.1 * ones(3); 0.001 * ones(3); 0.1 * ones(4); 0.01 * ones(3)] for i = 1:Nb]...)

ot1 = (x, u, w) -> transpose(x - zM) * Diagonal(Δt * qt) * (x - zM) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
ot2 = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT) + transpose(u) * Diagonal(Δt * 0.01 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(Δt * qt) * (x - zT)

ct1 = Cost(ot1, n, m, d)
ct2 = Cost(ot2, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]

# Constraints
function goal(x, u, w)
    Δ = x - zT
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
