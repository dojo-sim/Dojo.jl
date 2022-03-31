using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra 

# ## system
gravity = -9.81
timestep = 0.1
env = get_environment(:block, 
    representation=:maximal, 
    timestep=timestep,
    friction_coefficient=0.5,
    gravity=gravity);

# ## visualizer 
render(env.vis)

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z1 = [0.0; 0.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
zT = [0.0; 1.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # right goal
## zT = [0.0; 0.0; 0.25 + 1.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # up goal

# ## horizon
T = 11

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w, attitude_decompress=true),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w, attitude_decompress=true),
    n, n, m)
model = [dyn for t = 1:T-1];

# ## rollout
ū = [[0.0; 0.0; 0.0] for t = 1:T-1]
x̄ = rollout(model, z1, ū)
visualize(env, x̄);

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-2 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m)
cT = Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT];

# ## constraints
goal(x, u, w) = x - zT

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT];

# ## solver
s = IterativeLQR.solver(model, obj, cons, 
    opts=Options(
        linesearch=:armijo,
        α_min=1.0e-5,
        obj_tol=1.0e-3,
        grad_tol=1.0e-3,
        con_tol=0.005,
        max_iter=100,
        max_al_iter=10,
        ρ_init=1.0,
        ρ_scale=10.0,
        verbose=false))
IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄);

# ## solve
@time IterativeLQR.solve!(s);

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(s)
@show IterativeLQR.eval_obj(s.m_data.obj.costs, s.m_data.x, s.m_data.u, s.m_data.w)
@show s.s_data.iter[1]
@show norm(goal(s.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
z_vis = [[z_sol[1] for t = 1:10]..., z_sol..., [z_sol[end] for t = 1:10]...]
u_vis = [[u_sol[1] for t = 1:10]..., u_sol..., [u_sol[end] for t = 1:10]...]
vis, anim = visualize(env, z_vis)
vis, anim = visualize_force!(vis, anim, z_vis, u_vis) 

set_camera!(env.vis, 
    zoom=50.0, 
    cam_pos=[100.0, 0.0, 0.0]);
    
set_floor!(env.vis, 
    x=0.0, 
    y=4.0, 
    z=0.02, 
    color=RGBA(0.7, 0.7, 0.7, 1.0));


    
    
    