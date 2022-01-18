using IterativeLQR

# ## system
gravity = -9.81
dt = 0.1
env = make("cartpole", 
    mode=:min, 
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
z1 = max2min(env.mechanism, cartpole_nominal_max())
zT = max2min(env.mechanism, cartpole_goal_max())

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
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-3 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - zT
end

cont = Constraint()
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem 
prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
IterativeLQR.solve!(prob,
    max_al_iter=10,
    verbose=true)

x_sol, u_sol = get_trajectory(prob)
visualize(env, x_sol)