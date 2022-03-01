using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra 

# ## system
gravity=-9.81
dt = 0.1
env = get_environment(:block, 
    representation=:maximal, 
    timestep=dt,
    friction_coefficient=0.5,
    gravity=gravity)

# ## visualizer 
open(env.vis) 

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z1 = [0.0; 0.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
zT = [1.0; 0.0; 0.25; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # right goal
## zT = [0.0; 0.0; 0.25 + 1.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # up goal

# ## horizon
T = 11

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [[0.0; 0.0; 0.0] for t = 1:T-1]
x̄ = rollout(model, z1, ū)
visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-2 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m)
cT = Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
goal(x, u, w) = x - zT

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem 
prob = IterativeLQR.solver(model, obj, cons, 
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
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
@time IterativeLQR.solve!(prob)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
v, anim = visualize(env, [[z_sol[1] for t = 1:10]..., z_sol..., [z_sol[end] for t = 1:10]...])

set_camera!(env.vis, 
    zoom=50.0, 
    cam_pos=[100.0, 0.0, 0.0])
    
set_floor!(env.vis, 
    x=0.0, 
    y=4.0, 
    z=0.02, 
    color=RGBA(0.7, 0.7, 0.7, 1.0))

