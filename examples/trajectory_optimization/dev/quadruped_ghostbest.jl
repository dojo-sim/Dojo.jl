## Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
# using IterativeLQR

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


# System
gravity = -9.81
timestep = 0.05
friction_coefficient = 0.8
mech = get_mechanism(:quadruped, timestep=timestep, gravity=gravity, friction_coefficient = friction_coefficient, damper = 0.0, spring = 0.0)

# Dimensions
H = 20
n = minimal_dimension(mech)
m = control_dimension(mech)
d = 0
u_mask = [zeros(12,6) I(12)]

# Reference trajectory
xref = quadruped_trajectory(mech, r = 0.05, z = 0.25; Δx = -0.06, Δfront = 0.08, N = Int(H/2), Ncycles = 1)
zref = [minimal_to_maximal(mech, x) for x in xref]
storage = generate_storage(mech, zref)
visualize(mech, storage, vis = vis)
x1ref = xref[1]
z1ref = zref[1]

# Initial GHOST state
ϵ0 = 1e-2
mech = get_mechanism(:quadruped, timestep=timestep, gravity=gravity, friction_coefficient = friction_coefficient, damper = 10.0, spring = 300.0)
initialize!(mech, :quadruped)
set_state!(mech, z1ref)
set_spring_offset!(mech, x1ref)
@elapsed storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0, undercut = 1.5)
visualize(mech, storage, vis = vis)
ghost_altitude = get_minimal_state(mech)[3]
xghost = deepcopy(xref)
for i = 1:H
	xghost[i][3] = ghost_altitude
end

# Initial conditions, controls, disturbances
no_contact_mech = get_mechanism(:quadruped, timestep=timestep, gravity=gravity, friction_coefficient = friction_coefficient, damper = 5.0, spring = 0.0, contact = false)
mech = get_mechanism(:quadruped, timestep=timestep, gravity=gravity, friction_coefficient = friction_coefficient, damper = 5.0, spring = 0.0)
w = [zeros(d) for t = 1:H-1]
ughost = [inverse_control(no_contact_mech, xghost[i], xghost[i+1]) for i = 1:H-1]



# Model
ϵtol = 1e-5
function fd(y, x, u, w)
	z = step!(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)
	y .= copy(maximal_to_minimal(mech, z))
end
function fdx(fx, x, u, w)
	fx .= copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[1])
end
function fdu(fu, x, u, w)
	∇u = copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[2])
	fu .= ∇u
end

# Time
dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
model = [dyn for t = 1:H-1]

# # Rollout
# xrol = rollout(model, xghost[1], ughost, w)
# storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in x̄])
# visualize(mech, storage; vis = vis)
# plot(hcat(get_sdf(mech, storage)...))

xsol = deepcopy(xghost)
# xsol = deepcopy(xrol)
usol = deepcopy(ughost)
xabs = deepcopy(xref)

Xsol = [xsol]
Usol = [usol]
################################################################################
################################################################################


# Objective
qt = [0.3; 0.05; 0.002; 0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
ots = [(x, u, w) -> transpose(x - xabs[t]) * Diagonal(timestep * qt) * (x - xabs[t]) +
	transpose(u) * Diagonal(timestep * [0.001*ones(6); 0.01 * ones(12)]) * u for t = 1:H-1]|
oT = (x, u, w) -> transpose(x - xabs[end]) * Diagonal(timestep * qt) * (x - xabs[end])

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
cons = [[cont for t = 1:H-1]..., conT]

prob = problem_data(model, obj, cons)

initialize_controls!(prob, usol)
initialize_states!(prob, xsol)

# Solve
ghost_ilqr_solve!(prob,
    verbose = true,
	linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
	con_tol=1e-4,
    max_iter=100,
    max_al_iter=20,
    ρ_init=1.0,
    ρ_scale=2.0)

a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
a = 10

# for iter = 1:1
# 	# Model
# 	ϵtol = 1e-8
# 	undercut = Inf
# 	function fd(y, x, u, w)
# 		z = step!(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)
# 		y .= copy(maximal_to_minimal(mech, z))
# 	end
# 	function fdx(fx, x, u, w)
# 		fx .= copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)[1])
# 	end
# 	function fdu(fu, x, u, w)
# 		∇u = copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)[2])
# 		fu .= ∇u
# 	end
#
# 	# Objective
# 	qt = [0.3; 0.05; 0.002; 0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
# 	ots = [(x, u, w) -> transpose(x - xabs[t]) * Diagonal(timestep * qt) * (x - xabs[t]) +
# 		transpose(u) * Diagonal(timestep * [0.01*ones(6); 0.01 * ones(12)]) * u for t = 1:H-1]
# 	oT = (x, u, w) -> transpose(x - xabs[end]) * Diagonal(timestep * qt) * (x - xabs[end])
#
# 	cts = Cost.(ots, n, m, d)
# 	cT = Cost(oT, n, 0, 0)
# 	obj = [cts..., cT]
#
# 	# Constraints
# 	function goal(x, u, w)
# 	    Δ = x - xabs[end]
# 	    return Δ[[1:2; 7:9; 13:2:36] ]
# 	end
# 	function slack(x, u, w)
# 	    Δ = u[1:6]
# 	    return 0.1*Δ
# 	end
#
# 	# cont = Constraint()
# 	cont = Constraint(slack, n, m)
# 	conT = Constraint(goal, n, 0)
# 	cons = [[cont for t = 1:H-1]..., conT]
#
# 	prob = problem_data(model, obj, cons)
# 	initialize_controls!(prob, usol)
# 	initialize_states!(prob, xsol)
#
# 	# Solve
# 	constrained_ilqr_solve!(prob,
# 	    verbose = true,
# 		linesearch=:armijo,
# 	    α_min=1.0e-5,
# 	    obj_tol=1.0e-3,
# 	    grad_tol=1.0e-3,
# 	    max_iter=100,
# 	    max_al_iter=5,
# 	    ρ_init=1.0,
# 	    ρ_scale=3.0)
#
# 	xsol, usol = get_trajectory(prob)
# 	push!(Xsol, xsol)
# 	push!(Usol, usol)
# 	println("U violation: ",  scn(norm([norm(u[1:6], Inf) for u in usol], Inf)))
#
# 	println("Solution: ###################################")
# 	storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in xsol])
# 	visualize(mech, storage, vis = vis)
# 	sleep(5.0)
#
# 	println("Rollout: #####################################")
# 	xrol = rollout(model, x1, usol, w)
# 	storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in xrol])
# 	visualize(mech, storage; vis = vis)
# 	sleep(5.0)
# end

plot(hcat(Usol[end]...)')
plot(hcat(Xsol[end]...)')


plot(hcat([u[1:6] for u in Usol[end]]...)')
ustar = deepcopy(Usol[end])

visualize(mech, storage; vis = vis)


mech = get_mechanism(:quadruped, timestep=timestep, gravity=gravity, friction_coefficient = friction_coefficient, damper = 5.0, spring = 0.0)
initialize!(mech, :quadruped)
set_state!(mech, minimal_to_maximal(mech, xabs[1]))

function controller!(mechanism, k)
	@show k
	# set_input!(mechanism, u_mask' * u_mask * ustar[k])
	set_input!(mechanism, ū[k])
    return
end
@elapsed storage = simulate!(mech, 0.95, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in Xsol[2]])
visualize(mech, storage, vis = vis)




function plot_cone(s, Δs; plt = plot())
	θ = 0:0.03:2π
	Xc = cos.(θ)
	Yc = sin.(θ)
	plot!(Xc, Yc, legend = false, aspectratio = 1.0, xlims = (-2,2), ylims = (-2,2))
	scatter!(1/s[1]*s[2:2], 1/s[1]*s[3:3], markersize = 6.0)
	α = 0:0.05:1.0
	Sα = [s + αi*Δs for αi in α]
	Xα = [1/s[1]*s[2] for s in Sα]
	Yα = [1/s[1]*s[3] for s in Sα]
	plot!(Xα, Yα)
	display(plt)
	return nothing
end


s = [5.3e-2, 4.8e-2, -2.2e-2]
Δs = [-1.2e-1, -1.1e-1, 3.1e-2]

sc  =[ 5.257121057914604e-2,  4.791302263356767e-2, -2.163502818806898e-2]
Δsc =[-1.160143187010312e-1, -1.131270199417930e-1,  3.107149823590583e-2]
sp  =[ 5.257121057914604e-2,  4.791302263356767e-2, -2.163502818806898e-2]
Δsp =[-5.644876087925099e-2, -5.203584025551676e-2,  2.131470920823536e-2]

sp =[ 7.329911284149580e-6,  1.334490802446364e-6, -1.846287829175919e-7]
Δsp =[-7.346875843985361e-6, -1.333471393251204e-6,  1.795267308829321e-7]
sc =[ 7.329911284149580e-6,  1.334490802446364e-6, -1.846287829175919e-7]
Δsc =[-7.252246777012890e-6, -1.327343585041552e-6,  1.705722367500799e-7]
plt = plot()
plot_cone(sp, Δsp, plt = plt)
plot_cone(sc, Δsc, plt = plt)

using Plots

function rosenbrock(x::Vector)
  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
function barrier(s)
	return -0.5 * log(s[1]^2 - s[2:end]'*s[2:end])
end
function barrier_slice(x; s1 = 1.0)
	return -0.5 * log(max(s1^2 - x'*x, 1e-3))
end

default(size=(600,600), fc=:heat)
x, y = -1.5:0.1:1.5, -1.5:0.1:1.5
z = Surface((x,y)->barrier_slice([x,y]), x, y)
surface(x,y,z, linealpha = 0.3)



ϵtol = 1e-5
function fd(y, x, u, w)
	z = step!(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)
	y .= copy(maximal_to_minimal(mech, z))
end
function fdx(fx, x, u, w)
	fx .= copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[1])
end
function fdu(fu, x, u, w)
	∇u = copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = 1.5, verbose = false)[2])
	fu .= ∇u
end

println("Rollout: #####################################")
xrol = rollout(model[1:7], x1ref, Usol[end], w)
xrol = rollout(model, x1ref, Usol[end], w)
storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in xrol])
visualize(mech, storage; vis = vis)

plot(hcat(get_sdf(mech, storage)...))








function ghost_ilqr_solve!(prob::ProblemData;
    linesearch=:armijo,
    max_iter=10,
	max_al_iter=10,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
	con_tol=1.0e-3,
	con_norm_type=Inf,
	ρ_init=1.0,
	ρ_scale=10.0,
	ρ_max=1.0e8,
    verbose=false)

	ϵinit = 1e-2

	# reset duals
    for (t, λ) in enumerate(prob.m_data.obj.λ)
        fill!(λ, 0.0)
	end

	# initialize penalty
	for (t, ρ) in enumerate(prob.m_data.obj.ρ)
        fill!(ρ, ρ_init)
	end

	for i = 1:max_al_iter
		verbose && println("  al iter: $i")

		########################################################################
		########################################################################
		# update model
		ϵtol = ϵinit / 2^(i-1)
		undercut = Inf
		function fd(y, x, u, w)
			z = step!(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)
			y .= copy(maximal_to_minimal(mech, z))
		end
		function fdx(fx, x, u, w)
			fx .= copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)[1])
		end
		function fdu(fu, x, u, w)
			fu .= copy(get_minimal_gradients(mech, minimal_to_maximal(mech, x), u, ϵ = ϵtol, btol = ϵtol, undercut = undercut, verbose = false)[2])
		end
		# Time
		dyn = Dynamics(fd, fdx, fdu, n, n, m, d)
		model = [dyn for t = 1:H-1]
		prob.m_data.model = model
		########################################################################
		########################################################################

		# primal minimization
		ilqr_solve!(prob,
            linesearch=linesearch,
            α_min=α_min,
		    max_iter=max_iter,
            obj_tol=obj_tol,
		    grad_tol=grad_tol,
		    verbose=verbose)

		# update trajectories
		objective!(prob.s_data, prob.m_data, mode=:nominal)

        # constraint violation
		prob.s_data.c_max[1] <= con_tol && break

        # dual ascent
		augmented_lagrangian_update!(prob.m_data.obj,
			s=ρ_scale, max_penalty=ρ_max)

		########################################################################
		########################################################################
		xsol, usol = get_trajectory(prob)
		push!(Xsol, xsol)
		push!(Usol, usol)
		println("U violation: ",  scn(norm([norm(u[1:6], Inf) for u in usol], Inf)))

		println("Solution: ###################################")
		storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in xsol])
		visualize(mech, storage, vis = vis)
		sleep(5.0)

		println("Rollout: #####################################")
		xrol = rollout(model, x1ref, usol, w)
		storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in xrol])
		visualize(mech, storage; vis = vis)
		sleep(5.0)
		########################################################################
		########################################################################
	end

    return nothing
end
