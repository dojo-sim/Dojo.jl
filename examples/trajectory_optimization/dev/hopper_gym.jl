using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
dt = 0.05
gravity=-9.81
env = get_environment("hopper",
    timestep=dt,
    gravity=gravity);

# ## visualizer
open(env.vis)
render(env)

# ## dimensions
nx = env.num_states
nu = env.num_inputs
d = 0

# ## horizon 
T = 31
Tm = convert(Int, floor((T-1) / 2))

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    nx, nx, nu, d)
model = [dyn for t = 1:T-1] 

# ## initial conditions
x1 = [1.025; 0.0; 0.25; 0.0; 0.0; 0.0; -1.0; 0.0; 1.25; 0.0; -0.5; 0.0]
xM = [1.025 + 1.0; 0.0; 0.25; 0.0; 0.0; 0.0; -1.0; 0.0; 1.25; 0.0; -0.5; 0.0]
xT = [1.025 + 1.0; 0.0; 0.25; 0.0; 0.0; 0.0; -1.0; 0.0; 1.25; 0.0; -0.5; 0.0]

# # rollout
ū = [0.01 * randn(3) for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## object
obj1 = (x, u, w) -> 10.0 * transpose(x - xM) * Diagonal([1.0; 1.0; 1.0; 1.0e-3; 1.0e-3; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3]) * (x - xM) + transpose(u) * Diagonal(1.0e-1 * [1.0; 1.0; 0.1]) * u
obj2 = (x, u, w) -> 10.0 * transpose(x - xT) * Diagonal([1.0; 1.0; 1.0; 1.0e-3; 1.0e-3; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3]) * (x - xT) + transpose(u) * Diagonal(1.0e-1 * [1.0; 1.0; 0.1]) * u
objT = (x, u, w) -> 10.0 * transpose(x - xT) * Diagonal([1.0; 1.0; 1.0; 1.0e-3; 1.0e-3; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3; 1.0; 1.0e-3]) * (x - xT)

ct1 = IterativeLQR.Cost(obj1, nx, nu, 0)
ct2 = IterativeLQR.Cost(obj2, nx, nu, 0)
cT = IterativeLQR.Cost(objT, nx, 0, 0)
obj = [[ct1 for t = 1:Tm]..., [ct2 for t = 1:Tm]..., cT]

# ## constraints
goal(x, u, w) = (x - xT)[collect([1,2,3])]#,7,9,11])]

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# ## problem
prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
solve!(prob, 
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    con_tol=0.005,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0, 
    verbose=true)

# ## solution
x_sol, u_sol = get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## 
# using Plot
# plot(hcat(x_sol...)')
# plot(hcat(u_sol..., u_sol[end])', linetype=:steppost)

# ## visualize
x_vis = [[x̄[1] for t = 1:10]..., x̄..., [x̄[end] for t = 1:10]...]
x_vis = [[x_sol[1] for t = 1:10]..., x_sol..., [x_sol[end] for t = 1:10]...]

x_vis = [[1.025; 0.0; 0.25; 0.0; 0.0; 0.0; -1.0; 0.0; 1.25; 0.0; -0.5; 0.0]]
visualize(env, [[x_vis[1] for t = 1:10]..., x_vis..., [x_vis[end] for t = 1:10]...])

plot(hcat(u_sol...)')