using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# ## setup
using Dojo
using Plots
# using OSFLoader
using IterativeLQR
using LinearAlgebra
using DojoEnvironments

# ## visualizer
vis = Visualizer()
render(vis)
open(vis)

# using Pkg
# Pkg.develop(path="~/.julia/dev/Dojo.jl/DojoEnvironments")

# ## system
sphere_mass = 10.0
sphere_radius = 0.25
gravity = -9.81
timestep = 0.02
friction_coefficient = 0.05
collider_options = ColliderOptions(
	sliding_friction=friction_coefficient,
	impact_damper=3e5,
	impact_spring=3e4,
	)
################################################################################
# Simulation
################################################################################
mech = DojoEnvironments.get_mechanism(:nerf_sphere,
	nerf=:bunny,
	timestep=timestep,
	gravity=gravity,
    friction_coefficient=friction_coefficient,
	radius=sphere_radius,
	mass=sphere_mass,
    collider_options=collider_options)

# initial conditions
x_nerf = [1.12, -0.47, 0.05]
q_nerf = Quaternion(normalize([0.72, 0.692, -0.023, 0.033])..., false)
x_sphere = [2.0, -0.5, 0.0]
q_sphere = Quaternion(1.0, 0.0, 0.0, 0.0, false)

# initial state
z_nerf = [x_nerf; zeros(3); vector(q_nerf); zeros(3)]
z_sphere = [x_sphere; zeros(3); vector(q_sphere); zeros(3)]
z_initial = [z_nerf; z_sphere]
set_maximal_state!(mech, deepcopy(z_initial))

initialize!(mech, :nerf_sphere,
    nerf_position=x_nerf,
    nerf_orientation=q_nerf,
    sphere_position=x_sphere,
    sphere_orientation=q_sphere,
    )
z_initial = get_maximal_state(mech)

function control!(m, k; kp=1*1e0, kd=1*3e-1, xg_sphere=[-0,-0.5,0.5], xg_nerf=[-1.5, 0.40,0.369])
    sphere = get_body(m, :sphere)
    nerf = get_body(m, :bunny)
    x_sphere = current_position(sphere.state)
    v_sphere = current_velocity(sphere.state)[1]
    x_nerf = current_position(nerf.state)
    v_nerf = current_velocity(nerf.state)[1]
    u_sphere = kp * (xg_sphere - x_sphere) - kd * v_sphere
	@show u_sphere
    set_input!(m, [szeros(6); u_sphere; szeros(3)] / timestep)
    return nothing
end

storage = simulate!(mech, 3.20, control!, opts=SolverOptions(rtol=3e-4, btol=3e-4))
# final state
z_final = get_maximal_state(mech)
z_final[2] += 1.0
visualize(mech, generate_storage(mech, [z_initial]), vis=vis)
visualize(mech, generate_storage(mech, [z_final]), vis=vis)
visualize(mech, storage, vis=vis)
# current_position(mech.bodies[1].state)
# current_orientation(mech.bodies[1].state)

################################################################################
# Optimization
################################################################################
env = get_environment(:nerf_sphere,
    nerf=:bunny,
    vis=vis,
    representation=:minimal,
    friction_coefficient=friction_coefficient,
    timestep=timestep,
    gravity=gravity,
	radius=sphere_radius,
	mass=sphere_mass,
	collider_options=collider_options,
    infeasible_control=false, # it seems to work better without infeasible control because with nerf contact we can't make the contact `active at a distance` like we can with dojo-contact
    opts_step=SolverOptions(rtol=3e-3, btol=3e-3),
    opts_grad=SolverOptions(rtol=3e-3, btol=3e-3),
    )

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z_initial
z1 = deepcopy(z_initial)
x1 = maximal_to_minimal(env.mechanism, z_initial)
xT = maximal_to_minimal(env.mechanism, z_final)
x1 = maximal_to_minimal(mech, z_initial)
xT = maximal_to_minimal(mech, z_final)


# ## horizon
T = 60

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
# ū = [1.0 * [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0] for t = 1:T-1]
ū = [1.0 * [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0] for t = 1:T-1]

x̄ = IterativeLQR.rollout(model, x1, ū)
DojoEnvironments.visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT) +
    transpose(u) * Diagonal(2.0e-1 * ones(m)) * u
oT = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT)

ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]


# ## constraints
############################################################################
# ul = -1.0 * 1e-3*ones(nu_infeasible) #TODO this is relaxed
# uu = +1.0 * 1e-3*ones(nu_infeasible) #TODO this is relaxed
#
# function contt(x, u, w)
#     [
#         1e-1 * (ul - u[1:nu_infeasible]);
#         1e-1 * (u[1:nu_infeasible] - uu);
#     ]
# end

function goal(x, u, w)
	# x[[1,2,3,7,8,9,13,14,15]] - xT[[1,2,3,7,8,9,13,14,15]]
    1e-0 * (x[[1,2,3,13,14,15]] - xT[[1,2,3,13,14,15]])
    # x[[1]] - xT[[1]]
end

# con_policyt = IterativeLQR.Constraint(contt, n, m, indices_inequality=collect(1:2nu_infeasible))
con_policyt = IterativeLQR.Constraint()
con_policyT = IterativeLQR.Constraint(goal, n, 0)

cons = [[con_policyt for t = 1:T-1]..., con_policyT]

# ## solver
solver_options = IterativeLQR.Options(
    line_search=:armijo,
	# max_iterations=75,
    max_iterations=1,
	# max_dual_updates=10,
    max_dual_updates=1,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    scaling_penalty=10.0,
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
scale_vis = 0.25
offset = [-0.75, 0.0, 0.0]
for body in env.mechanism.bodies
	body.shape.scale *= scale_vis
end
sol_storage = generate_storage(env.mechanism, [minimal_to_maximal(env.mechanism, z) for z in z_sol])
for x in sol_storage.x
	x .*= scale_vis
	x .+= fill(offset, 60)
end
vis, anim = visualize(env.mechanism, sol_storage, vis=vis, color=RGBA(0.9, 0.9, 0.9, 1.0))



################################################################################
# Export
################################################################################
# vis = Visualizer()
# render(vis)

# # initial guess
# sol_storage = generate_storage(mech, [minimal_to_maximal(mech, x) for x in x̄])
# for x in sol_storage.x
# 	x .*= scale_vis
# 	x .+= fill(offset, 60)
# end


padded_storage = generate_storage(mech, [
	fill(get_maximal_state(sol_storage, 1), 30);
	get_maximal_state(sol_storage);
	fill(get_maximal_state(sol_storage, T), 30)])
vis, anim = visualize(env.mechanism, padded_storage, vis=vis, color=RGBA(0.9, 0.9, 0.9, 1.0))

panda_mech = get_mechanism(:panda, damper=100.0, model_type=:end_effector, contact=true,
	joint_limits=[[-10.0, -1.7628, -2.8973, -0.0698, -2.8973, -3.7525, -2.8973, -0.00],
				  [ 10.0,  1.7628,  2.8973,  3.0718,  2.8973,  0.0175,  2.8973,  0.04]],)

initialize!(panda_mech, :panda, joint_angles=[π, -0.8, 0.0, 1.6, 0.0, -3.2, 0.0, 0.0, 0.0])
storage = simulate!(panda_mech, 0.01)
visualize(panda_mech, storage, vis=vis, name=:panda)
q_end_effector = current_orientation(get_body(panda_mech, 12).state)

visual_storage = deepcopy(padded_storage)

z_panda = panda_inverse_kinematics_trajectory(panda_mech,
	visual_storage.x[2][1:T+60],
	[q_end_effector for i=1:T+60])

panda_storage = generate_storage(panda_mech, z_panda)
vis, anim = visualize(panda_mech, panda_storage, vis=vis, animation=anim, name=:panda,
	show_contact=true)

function set_target!(vis::Visualizer)
	cyl = MeshCat.Cylinder(Point(0,0,0.0), Point(0,0,0.002), 0.1)
	setobject!(vis[:sphere_target], cyl, MeshPhongMaterial(color=RGBA(0.9,0.9,0.9,1.0)))
	setobject!(vis[:nerf_target], cyl, MeshPhongMaterial(color=RGBA(0.2,0.2,0.2,1.0)))
	settransform!(vis[:sphere_target], Translation([xT[1:2]; 0]*scale_vis + offset))
	settransform!(vis[:nerf_target], Translation([xT[13:14]; 0]*scale_vis + offset))
	return nothing
end
set_target!(vis)

convert_frames_to_video_and_gif("sphere_bunny_trajopt_al0_padded")
