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

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


# System
gravity = -9.81
Δt = 0.05
cf = 0.8
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = cf, damper = 0.0, spring = 0.0)

# Dimensions
T = 20
n = minCoordDim(mech)
m = controldim(mech)
d = 0
u_mask = [zeros(12,6) I(12)]

# Reference trajectory
xref = quadruped_trajectory(mech, r = 0.05, z = 0.25; Δx = -0.06, Δfront = 0.08, N = Int(T/2), Ncycles = 1)
zref = [min2max(mech, x) for x in xref]
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)
x1 = xref[1]
z1 = zref[1]

# Gravity compensation
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = cf, damper = 10.0, spring = 300.0)
initialize!(mech, :quadruped)
setState!(mech, z1)
setSpringOffset!(mech, x1)
@elapsed storage = simulate!(mech, 4.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
ugc = gravity_compensation(mech)



# Initial conditions, controls, disturbances
mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = cf, damper = 5.0, spring = 0.0)
u_control = u_mask * ugc
ū = [0.9*u_mask'*u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]

# Model
ϵtol = 3e-2
function fd(y, x, u, w)
	z = simon_step!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)
	y .= copy(max2min(mech, z))
end
function fdx(fx, x, u, w)
	fx .= copy(getMinGradients!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[1])
end
function fdu(fu, x, u, w)
	∇u = copy(getMinGradients!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[2])
	fu .= ∇u
end

# Time
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:T-1]

# Rollout
x̄ = rollout(model, x1, ū, w)
storage = generate_storage(mech, [min2max(mech, x) for x in x̄])
visualize(mech, storage; vis = vis)



xsol = deepcopy(x̄)
usol = deepcopy(ū)
xabs = deepcopy(xref)

Xsol = [xsol]
Usol = [usol]
################################################################################
################################################################################


for iter = 1:1
	# Model
	ϵtol = 1e-5
	function fd(y, x, u, w)
		z = simon_step!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)
		y .= copy(max2min(mech, z))
	end
	function fdx(fx, x, u, w)
		fx .= copy(getMinGradients!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[1])
	end
	function fdu(fu, x, u, w)
		∇u = copy(getMinGradients!(mech, min2max(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[2])
		fu .= ∇u
	end

	# Objective
	qt = [0.3; 0.05; 0.01; 0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
	ots = [(x, u, w) -> transpose(x - xabs[t]) * Diagonal(Δt * qt) * (x - xabs[t]) +
		transpose(u) * Diagonal(Δt * [0.01*ones(6); 0.01 * ones(12)]) * u for t = 1:T-1]
	oT = (x, u, w) -> transpose(x - xabs[end]) * Diagonal(Δt * qt) * (x - xabs[end])

	cts = Cost.(ots, n, m, d)
	cT = Cost(oT, n, 0, 0)
	obj = [cts..., cT]

	# Constraints
	function goal(x, u, w)
	    Δ = x - xabs[end]
	    return Δ[[1:2; 7:9; 13:2:36] ]
	end
	function slack(x, u, w)
	    Δ = u[1:6]
	    return 0.1*Δ
	end

	# cont = Constraint()
	cont = Constraint(slack, n, m)
	conT = Constraint(goal, n, 0)
	cons = [[cont for t = 1:T-1]..., conT]

	prob = problem_data(model, obj, cons)
	initialize_controls!(prob, usol)
	initialize_states!(prob, xsol)

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

	xsol, usol = get_trajectory(prob)
	push!(Xsol, xsol)
	push!(Usol, usol)
	println("U violation: ",  scn(norm([norm(u[1:6], Inf) for u in usol], Inf)))

	println("Solution: ###################################")
	storage = generate_storage(mech, [min2max(mech, x) for x in xsol])
	visualize(mech, storage, vis = vis)
	sleep(5.0)

	println("Rollout: #####################################")
	xrol = rollout(model, x1, usol, w)
	storage = generate_storage(mech, [min2max(mech, x) for x in xrol])
	visualize(mech, storage; vis = vis)
	sleep(5.0)
end

plot(hcat(Usol[end]...)')
plot(hcat(Xsol[end]...)')

ustar = deepcopy(Usol[end])


mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = cf, damper = 5.0, spring = 0.0)
initialize!(mech, :quadruped)
setState!(mech, min2max(mech, xabs[1]))

function controller!(mechanism, k)
	@show k
	setControl!(mechanism, u_mask' * u_mask * ustar[k])
    return
end
@elapsed storage = simulate!(mech, 0.95, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
