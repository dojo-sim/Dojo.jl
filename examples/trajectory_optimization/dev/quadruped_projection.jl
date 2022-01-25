# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

using IterativeLQR

# System
gravity = -9.81
Δt = 0.05
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = 0.8, damper = 10.0, spring = 0.0)
initialize!(mech, :quadruped, tran = [0,0,0.], v = [0.5,0,0.])
# x0 = getMinState(mech)
# x0[35] = 0.4
# set_state!(mech, min2max(mech, x0))
@elapsed storage = simulate!(mech, 0.05, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

n = minimal_dimension(mech)
m = 12 + n
d = 0
T = 18

xref = quadruped_trajectory(mech, r = 0.08, z = 0.27; Δt = Δt, Δx = -0.04, Δfront = 0.10, N = Int(T/2), Ncycles = 1)
zref = [min2max(mech, x) for x in xref]
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)
x1 = xref[1]
z1 = zref[1]
# visualizeMaxCoord(mech, zref[1], vis)

function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = control_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for eqc in mechanism.eqconstraints
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

mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = 0.8, damper = 100.0, spring = 200.0)
initialize!(mech, :quadruped)
set_state!(mech, z1)
setSpringOffset!(mech, x1)
@elapsed storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
ugc = gravity_compensation(mech)

mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = 0.8, damper = 2.0, spring = 0.0)
u_control = ugc[6 .+ (1:12)]
u_mask = [zeros(12,6) I(12)]

z = [copy(z1)]
for t = 1:5
    znext = step!(mech, z[end], u_mask'*u_control)
    push!(z, znext)
end
storage = generate_storage(mech, z)
visualize(mech, storage, vis = vis)


# Model
function fd(y, x, u, w)
    u_control = u[1:12]
    s = u[12 .+ (1:n)]
	z = step!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)
	yvio = copy(max2min(mech, z)) + s
	y .= projectQuadruped!(mech, yvio)
end

function fdx(fx, x, u, w)
	u_control = u[1:12]
    s = u[12 .+ (1:n)]
	fx .= copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[1])
end

function fdu(fu, x, u, w)
	u_control = u[1:12]
    s = u[12 .+ (1:n)]
	∇u = copy(getMinGradients!(mech, min2max(mech, x), u_mask'*u_control, ϵ = 3e-4, btol = 3e-4, undercut = 1.5, verbose = false)[2])
	fu .= [∇u * u_mask' I(n)]
end


# Time
h = mech.Δt
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]


# Initial conditions, controls, disturbances
# ū = [[u_control; zeros(n)] for t = 1:T-1]
# ū = [[u_control; xref[t+1] - xref[t]] for t = 1:T-1]
# ū = [[zeros(12); xref[t+1] - xref[t]] for t = 1:T-1]
ū = [[u_control; xref[t+1] - xref[t]] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Rollout
# x̄ = rollout(model, x1, ū, w)
x̄ = deepcopy(xref)
storage = generate_storage(mech, [min2max(mech, x) for x in x̄])
visualize(mech, storage; vis = vis)

# Objective
qt = [0.3; 0.1; 0.1; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(Δt * qt) * (x - xref[t]) +
	transpose(u) * Diagonal(Δt * [0.01*ones(12); 10*ones(n)]) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(Δt * qt) * (x - xref[end])

cts = Cost.(ots, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [cts..., cT]

# Constraints
function goal(x, u, w)
    Δ = x - xref[end]
    return Δ[collect(1:3)]
end
function slack(x, u, w)
    Δ = u[12 .+ (1:n)]
    return 0.01*Δ
end

cont = Constraint()
# cont = Constraint(slack, n, m)
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
    ρ_scale=3.0)

x_sol, u_sol = get_trajectory(prob)
storage = generate_storage(mech, [min2max(mech, x) for x in x_sol])
visualize(mech, storage, vis = vis)

norm([norm(u[12 .+ (1:n)], Inf) for u in u_sol], Inf)



function projectQuadruped!(mechanism::Mechanism{T}, x::AbstractVector{T}) where {T}
	for leg in [:RL, :RR, :FL, :FR]
		x = projectQuadrupedLeg!(mechanism, x; leg = leg)
	end
	xp = getMinState(mechanism)
	return xp
end

function projectQuadrupedLeg!(mechanism::Mechanism{T}, x::AbstractVector{T}; leg::Symbol = :FR) where {T}
	xp = x
	z = min2max(mechanism, x)
	set_state!(mech, z)

	# starting point of the local search
	θhip = minimal_coordinates(mech, get_joint_constraint(mech, String(leg)*"_thigh_joint"))
	θknee = minimal_coordinates(mech, get_joint_constraint(mech, String(leg)*"_calf_joint"))
	θ = [θhip; θknee]
	for k = 1:10
		s = sdfquadruped(mechanism, θ; leg = leg)
		(s > -1e-10) && continue
		∇ = FiniteDiff.finite_difference_jacobian(θ -> sdfquadruped(mechanism, θ; leg = leg), θ)
		θ -= ∇ \ [s]
	end
	xp = getMinState(mechanism)
	return xp
end

function sdfquadruped(mechanism::Mechanism{T}, θ::AbstractVector{T}; leg::Symbol = :FR) where {T}
	set_position(mechanism, get_joint_constraint(mechanism, String(leg)*"_thigh_joint"), [θ[1]])
	set_position(mechanism, get_joint_constraint(mechanism, String(leg)*"_calf_joint"), [θ[2]])

	foot = get_body(mechanism, String(leg)*"_calf")
	ineqcs = collect(mechanism.ineqconstraints)
	ineqc = ineqcs[findfirst(x -> x.parentid == foot.id, ineqcs)]
	p = contact_location(ineqc, foot)
	return p[3]
end



# zgood = deepcopy(z0)
# z0 = deepcopy(zgood)
# set_state!(mech, z0)
# z0 = getState(mech)
x0 = max2min(mech, z0)
x1 = deepcopy(x0)
x1[3] -= 0.1
z1 = min2max(mech, x1)
# visualizeMaxCoord(mech, z0, vis)
visualizeMaxCoord(mech, z1, vis)
xp1 = projectQuadruped!(mech, x1)
zp1 = min2max(mech, xp1)
set_state!(mech, zp1)
visualizeMaxCoord(mech, zp1, vis)

contact_location(mech)
