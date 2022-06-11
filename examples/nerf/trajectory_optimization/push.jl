using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# ## setup
using Dojo
using Plots
using OSFLoader
using IterativeLQR
using LinearAlgebra
using DojoEnvironments

# ## visualizer
vis = Visualizer()
render(vis)

# using Pkg
# Pkg.develop(path="~/.julia/dev/Dojo.jl/DojoEnvironments")

# ## system
gravity = -9.81
timestep = 0.02
friction_coefficient = 0.30
################################################################################
# Simulation
################################################################################
mech = DojoEnvironments.get_mechanism(:nerf_sphere, nerf=:bluesoap, timestep=timestep, gravity=gravity,
# mech = get_nerf_sphere(nerf=:bunny, timestep=timestep, gravity=gravity,
    friction_coefficient=friction_coefficient,
    collider_options=ColliderOptions(sliding_friction=friction_coefficient))

mech.contacts[end].model.collision.collider.options.impact_damper = 10000.0
mech.contacts[end].model.collision.collider.options.impact_spring = 1000.0

# initial conditions
x_bluesoap = [0.269, -0.40, 0.369]
# x_bluesoap = [0., -0.5, 0.45]
q_bluesoap = Quaternion(-0.247, 0.715, 0.618, -0.211, false)
# q_bluesoap = Quaternion(normalize([1, 0, 0, 0.5])..., false)
x_sphere = [2.0, -0.5, 0.5]
q_sphere = Quaternion(1.0, 0.0, 0.0, 0.0, false)

# initial state
z_bluesoap = [x_bluesoap; zeros(3); vector(q_bluesoap); zeros(3)]
z_sphere = [x_sphere; zeros(3); vector(q_sphere); zeros(3)]
z_initial = [z_bluesoap; z_sphere]
set_maximal_state!(mech, deepcopy(z_initial))


initialize!(mech, :nerf_sphere,
    nerf_position=x_bluesoap - [0,0,0.5],
    nerf_orientation=q_bluesoap,
    sphere_position=x_sphere - [0,0,0.5],
    sphere_orientation=q_sphere,
    )

function ctrl!(m, k; kp=1e-0, kv=3e-1, xg_sphere=[-0,-0.5,0.5], xg_nerf=[-1.5, 0.40,0.369])
    sphere = get_body(m, :sphere)
    nerf = get_body(m, :bluesoap)
    x_sphere = current_position(sphere.state)
    v_sphere = current_velocity(sphere.state)[1]
    x_nerf = current_position(nerf.state)
    v_nerf = current_velocity(nerf.state)[1]
    u_sphere = (xg_sphere - x_sphere) * kp - kv * v_sphere
	# u_nerf = szeros(3) # (xg_nerf - x_nerf) * kp - kv * v_nerf
    u_nerf = (xg_nerf - x_nerf) * kp - kv * v_nerf
    # set_input!(m, [u_nerf; szeros(3); u_sphere; szeros(3)])
    set_input!(m, [u_nerf; szeros(3); u_sphere; szeros(3)])
    return nothing
end

storage = simulate!(mech, 0.60, ctrl!, opts=SolverOptions(rtol=3e-4, btol=3e-4))
# final state
z_final = get_maximal_state(mech)
# z_final = deepcopy(z_initial)
# z_final[[1]] .+= 1.0
# z_final[[1,2,3]] .+= 1.0
# z_final[[14,15,16]] .+= 1.0
visualize(mech, storage, vis=vis)

visualize(mech, generate_storage(mech, [z_initial]), vis=vis)
visualize(mech, generate_storage(mech, [z_final]), vis=vis)


################################################################################
# Optimization
################################################################################
env = get_environment(:nerf_sphere,
    nerf=:bluesoap,
    vis=vis,
    representation=:minimal,
    # representation=:maximal,
    friction_coefficient=friction_coefficient,
    timestep=timestep,
    gravity=gravity,
    infeasible_control=true,
    opts_step=SolverOptions(rtol=3e-3, btol=3e-3),
    opts_grad=SolverOptions(rtol=3e-3, btol=3e-3),
    )

env.mechanism.contacts[end].model.collision.collider.options.impact_damper = 10000.0
env.mechanism.contacts[end].model.collision.collider.options.impact_spring = 1000.0


# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z_initial
z1 = deepcopy(z_initial)
x1 = maximal_to_minimal(env.mechanism, z_initial)
xT = maximal_to_minimal(env.mechanism, z_final)


# ## horizon
T = 30

function finite_difference_dynamics_jacobian_state(dx, env, x, u, w)
	function local_dynamics(x)
		y = zeros(24)
		dynamics(y, env, x, u, w)
		return y
	end
	dx .= FiniteDiff.finite_difference_jacobian(x -> local_dynamics(x), x)
	return nothing
end

function finite_difference_dynamics_jacobian_input(du, env, x, u, w)
	function local_dynamics(u)
		y = zeros(24)
		dynamics(y, env, x, u, w)
		return y
	end
	du .= FiniteDiff.finite_difference_jacobian(u -> local_dynamics(u), u)
	return nothing
end

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)
# dyn = IterativeLQR.Dynamics(
#     (y, x, u, w) -> dynamics(y, env, x, u, w),
#     (dx, x, u, w) -> finite_difference_dynamics_jacobian_state(dx, env, x, u, w),
#     (du, x, u, w) -> finite_difference_dynamics_jacobian_input(du, env, x, u, w),
#     n, n, m)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [1.0 * [-0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū)
DojoEnvironments.visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT) +
    transpose(u) * Diagonal(1.0e-1 * ones(m)) * u
oT = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT)

ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
function goal(x, u, w)
	# x[[1,2,3,7,8,9,13,14,15]] - xT[[1,2,3,7,8,9,13,14,15]]
    1e-0 * (x[[1,2,3,13,14,15]] - xT[[1,2,3,13,14,15]])
    # x[[1]] - xT[[1]]
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## solver
solver_options = IterativeLQR.Options(
    line_search=:armijo,
    max_iterations=75,
    max_dual_updates=16,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    scaling_penalty=3.0,
    max_penalty=1e4,
    verbose=true)
s = IterativeLQR.Solver(model, obj, cons, options=solver_options)
IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄)

# ## callback
function local_callback(solver; )
    u_sol = s.problem.actions
    x_sol = s.problem.states
    x_rollout = IterativeLQR.rollout(solver.problem.model.dynamics, x_sol[1], u_sol)
    DojoEnvironments.visualize(env, x_rollout)
    return nothing
end


# ## solve
@time IterativeLQR.constrained_ilqr_solve!(s,
	augmented_lagrangian_callback! = local_callback)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(s)

# ## visualize
DojoEnvironments.visualize(env, z_sol)


# convert_frames_to_video_and_gif("bluesoap_push_high_friction")
# obj = MeshFileGeometry(joinpath("/home/simon/Downloads/UM2_logo_7.obj"))
# setobject!(vis[:logo], obj)
# settransform!(vis[:logo], LinearMap(0.01*I))

# Dojo.convert_video_to_gif("/home/simon/Downloads/video_trial.mp4",
#     "/home/simon/Downloads/video_trial.gif", framerate=15, width=400)




################################################################################
# gradients coming from nerf object interactions are wrong
################################################################################





mech.bodies







nu = input_dimension(mech)
u1 = 1 * ones(nu)

initialize!(mech, :nerf_sphere,
    nerf_position=[0,0,1],
    sphere_position=[0,0,3],
    )
z1 = get_maximal_state(mech)
fx, fu = get_maximal_gradients!(env.mechanism, z1, u1, opts=env.opts_grad)
plot(Gray.(1e3abs.(fu)))
plot(Gray.(1e3abs.(fx)))
# fx, fu = get_minimal_gradients!(env.mechanism, z1, u1, opts=env.opts_grad)
# plot(Gray.(1e0abs.(fu)))
# plot(Gray.(1e0abs.(fx)))

function get_maximal_gradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	nu = input_dimension(mechanism)

	for entry in mechanism.data_matrix.nzval # reset matrix
		entry.value .= 0.0
	end
	jacobian_data!(mechanism.data_matrix, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = data_dim.(nodes)
	index_row = [1+sum(dimrow[1:i-1]):sum(dimrow[1:i]) for i in 1:length(dimrow)]
	index_col = [1+sum(dimcol[1:i-1]):sum(dimcol[1:i]) for i in 1:length(dimcol)]

	index_state = [index_col[body.id][[14:16; 8:10; 17:19; 11:13]] for body in mechanism.bodies] # ∂ x2 v15 q2 ϕ15
	index_control = [index_col[joint.id][1:input_dimension(joint)] for joint in mechanism.joints] # ∂ u

	datamat = full_matrix(mechanism.data_matrix, dimrow, dimcol)
	solmat = full_matrix(mechanism.system)

	# data Jacobian
	data_jacobian = solmat \ datamat #TODO: use pre-factorization

	# Jacobian
	jacobian_state = zeros(12Nb,12Nb)
	jacobian_control = zeros(12Nb,nu)
	for (i, body) in enumerate(mechanism.bodies)
		id = body.id
		# Fill in gradients of v25, ϕ25
		jacobian_state[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_state...)]
		jacobian_control[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_control...)]

		# Fill in gradients of x3, q3
		x2 = body.state.x2
		q2 = body.state.q2
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		q3 = next_orientation(q2, ϕ25, timestep)
		jacobian_state[12*(i-1) .+ (1:3), :] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (1:3), 12*(i-1) .+ (1:3)] += linear_integrator_jacobian_position(x2, v25, timestep)
		jacobian_state[12*(i-1) .+ (7:9), :] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (7:9), 12*(i-1) .+ (7:9)] += LVᵀmat(q3)' * rotational_integrator_jacobian_orientation(q2, ϕ25, timestep, attjac=true)

		jacobian_control[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_control...)]
		jacobian_control[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_control...)]
	end

	return jacobian_state, jacobian_control
end


initialize!(mech, :nerf_sphere,
    nerf_position=[0,0,1],
    sphere_position=[0,0,2],
    )

# vector ordering
mech.bodies[1] # nerf
mech.bodies[2] # sphere
# id ordering
mech.bodies[1].id # nerf
mech.bodies[2].id # sphere
# maximal state ordering
get_maximal_state(mech) # nerf, sphere
# minimal state ordering
get_minimal_state(mech) # nerf, sphere
# root to leaves
root_to_leaves_ordering(mech) # nerf, sphere
# inputs
mech.joints[1]
mech.joints[2]

jacobian_data!(mech.data_matrix, mech)

nodes = [mech.joints; mech.bodies; mech.contacts]
dimrow = length.(nodes)
dimcol = data_dim.(nodes)
index_row = [1+sum(dimrow[1:i-1]):sum(dimrow[1:i]) for i in 1:length(dimrow)]
index_col = [1+sum(dimcol[1:i-1]):sum(dimcol[1:i]) for i in 1:length(dimcol)]

index_state = [index_col[body.id][[14:16; 8:10; 17:19; 11:13]] for body in mech.bodies] # ∂ x2 v15 q2 ϕ15
index_control = [index_col[joint.id][1:input_dimension(joint)] for joint in mech.joints] # ∂ u

datamat = full_matrix(mech.data_matrix, dimrow, dimcol)
solmat = full_matrix(mech.system)

plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(solmat)))

# data Jacobian
data_jacobian = solmat \ datamat
plot(Gray.(abs.(data_jacobian)))

data_jacobian[index_row[3], vcat(index_control...)]
data_jacobian[index_row[4], vcat(index_control...)]
index_row[3]
index_row[4]
index_control


# z1 = minimal_to_maximal(env.mechanism, x1)
# zT = minimal_to_maximal(env.mechanism, xT)
#
# @show norm((z1 - z_initial)[1:13])
# @show norm((zT - z_final)[1:13])
#
# @show norm((z1 - z_initial)[4:7])
# @show norm((zT - z_final)[4:7])
#
# @show norm((z1 - z_initial)[8:10])
# @show norm((zT - z_final)[8:10])
#
# @show norm((z1 - z_initial)[11:13])
# @show norm((zT - z_final)[11:13])
#
# @show norm((z1 - z_initial)[14:end])
# @show norm((zT - z_final)[14:end])

z1 = deepcopy(z_initial)
x1 = maximal_to_minimal(env.mechanism, z_initial)
nu = input_dimension(env.mechanism)
nx = minimal_dimension(env.mechanism)
w1 = zeros(0)
u1 = zeros(nu)
u1[1] += 100.0
u1[2] += 100.0
u1[3] += 100.0
u1[4] += 100.0
u1[5] += 100.0
u1[6] += 100.0

u1[7] += 100.0
u1[8] += 100.0
u1[9] += 100.0
u1[10] += 100.0
u1[11] += 100.0
u1[12] += 100.0
y1 = zeros(nx)
dynamics(y1, env, x1, u1, w1)

nz = maximal_dimension(env.mechanism)
y1 = zeros(nz)
dynamics(y1, env, z1, u1, w1)
# z2 = minimal_to_maximal(env.mechanism, y1)

# dx1 = zeros(nx,nx)
du1 = zeros(nz,nu)
dynamics_jacobian_state(dx1, env, x1, u1, w1)
dynamics_jacobian_input(du1, env, z1, u1, w1)
plot(Gray.(1e3abs.(dx1)))
plot(Gray.(1e3abs.(du1)))

norm(du1)


vis, anim = visualize(mech, generate_storage(mech, [z1]), vis=vis, name=:initial)
vis, anim = visualize(mech, generate_storage(mech, [z2]), vis=vis, animation=anim, name=:final)
z_initial





impulse_map_jacobian_configuration
