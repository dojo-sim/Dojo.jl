using Dojo
using IterativeLQR
using LinearAlgebra 

# ## system
gravity = -9.81
dt = 0.1
env = make("cartpole", 
    mode=:max, 
    dt=dt,
    g=gravity);

mujoco_inertia!(env.mech)

# ## visualizer 
open(env.vis) 

# ## dimensions
n = env.nx
m = env.nu
d = 0

# ## states
z1 = cartpole_nominal_max()
zT = cartpole_goal_max()

# ## horizon
T = 26

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [t < 5 ? 1.0 * rand(m) : (t < 10 ? -1.0 * rand(m) : zeros(m)) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = rollout(model, z1, ū, w)
visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(dt * ones(n)) * (x - zT) + transpose(u) * Diagonal(dt * [0.1]) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(100.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - zT
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem 
prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
IterativeLQR.solve!(prob,
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0,
    verbose=true)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(prob)
visualize(env, [[z_sol[1] for t = 1:10]..., z_sol..., [z_sol[end] for t = 1:10]...])

# ## ghost
ghost(env, z_sol)

