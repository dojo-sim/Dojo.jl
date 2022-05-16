using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()

# ## setup
using Dojo
using Plots
using OSFLoader
using IterativeLQR
using LinearAlgebra

# ## visualizer
vis = Visualizer()
open(vis)

# ## system
gravity = -9.81
timestep = 0.01
friction_coefficient = 0.01
################################################################################
# Simulation
################################################################################
mech = get_nerf_sphere(nerf=:bluesoap, timestep=timestep, gravity=gravity,
    friction_coefficient=friction_coefficient,
    collider_options=ColliderOptions(sliding_friction=friction_coefficient))

# initial conditions
x_bluesoap = [0.269, -0.40, 0.369]
q_bluesoap = Quaternion(-0.247, 0.715, 0.618, -0.211, false)
x_sphere = [2.0, -0.5, 0.5]
q_sphere = Quaternion(1.0, 0.0, 0.0, 0.0, false)

# initial state
z_bluesoap = [x_bluesoap; zeros(3); vector(q_bluesoap); zeros(3)]
z_sphere = [x_sphere; zeros(3); vector(q_sphere); zeros(3)]
z_initial = [z_bluesoap; z_sphere]
set_maximal_state!(mech, z_initial)
x_initial = maximal_to_minimal(mech, z_initial)
set_minimal_state!(mech, x_initial)


# initialize!(mech, :nerf_sphere,
#     nerf_position=x_bluesoap - [0,0,0.5],
#     nerf_orientation=q_bluesoap,
#     sphere_position=x_sphere - [0,0,0.5],
#     sphere_orientation=q_sphere,
#     )

function ctrl!(m, k; kp=1e-0, kv=3e-1, xg_sphere=[-2,-0.5,0.5], xg_nerf=[-2,-0.40,0.369])
    sphere = get_body(m, :sphere)
    # nerf = get_body(m, :bluesoap)
    x_sphere = current_position(sphere.state)
    v_sphere = current_velocity(sphere.state)[1]
    # x_nerf = current_position(nerf.state)
    # v_nerf = current_velocity(nerf.state)[1]
    u_sphere = (xg_sphere - x_sphere) * kp - kv * v_sphere
    u_nerf = szeros(3) # (xg_nerf - x_nerf) * kp - kv * v_nerf
    set_input!(m, [u_nerf; szeros(3); u_sphere; szeros(3)])
    return nothing
end

storage = simulate!(mech, 10.0, ctrl!, opts=SolverOptions(rtol=3e-4, btol=3e-4))
# final state
z_final = get_maximal_state(mech)
z_final = deepcopy(z_initial)
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
    friction_coefficient=friction_coefficient,
    timestep=timestep,
    gravity=gravity);

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
x1 = maximal_to_minimal(env.mechanism, z_initial)
xT = maximal_to_minimal(env.mechanism, z_final)

# ## horizon
T = 200

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [1.0 * [-0.2, 0.0, 0.0, 0.0, 0.0, 0.0] for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū)
visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT) +
    transpose(u) * Diagonal(1.0e-3 * ones(m)) * u
oT = (x, u, w) -> transpose(x - xT) * Diagonal(1.0e-1 * ones(n)) * (x - xT)

ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
function goal(x, u, w)
    x[[1,2,3,7,8,9,13,14,15]] - xT[[1,2,3,7,8,9,13,14,15]]
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## solver
solver_options = IterativeLQR.Options(
    line_search=:armijo,
    max_iterations=75,
    max_dual_updates=8,
    objective_tolerance=1e-3,
    lagrangian_gradient_tolerance=1e-3,
    constraint_tolerance=1e-3,
    scaling_penalty=10.0,
    max_penalty=1e7,
    verbose=true)
s = IterativeLQR.Solver(model, obj, cons, options=solver_options)
IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄)

# ## callback
function local_callback(solver; )
    u_sol = s.problem.actions
    x_sol = s.problem.states
    x_rollout = IterativeLQR.rollout(solver.problem.model.dynamics, x_sol[1], u_sol)
    visualize(env, x_rollout, vis=env.vis)
    return nothing
end

# ## solve
@time IterativeLQR.constrained_ilqr_solve!(s, augmented_lagrangian_callback! = local_callback)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(s)

# ## visualize
visualize(env, z_sol)


# convert_frames_to_video_and_gif("bluesoap_push_high_friction")
# obj = MeshFileGeometry(joinpath("/home/simon/Downloads/UM2_logo_7.obj"))
# setobject!(vis[:logo], obj)
# settransform!(vis[:logo], LinearMap(0.01*I))
