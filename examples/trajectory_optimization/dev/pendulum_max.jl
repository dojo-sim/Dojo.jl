# Open visualizer
vis = Visualizer()
open(vis)

using IterativeLQR

# system
dt = 0.1
gravity = -9.81
env = make("Pendulum", 
    mode=:max, 
    dt=dt,
    g=gravity,
    vis=vis);

## state space
n = env.nx
m = env.nu
d = 0

## initial state
z1 = pendulum_nominal_max()

## goal state 
xT = [0.0; 0.0; 0.5]
vT = [0.0; 0.0; 0.0]
qT = [0.0; 1.0; 0.0; 0.0]
ωT = [0.0; 0.0; 0.0]
zT = [xT; vT; qT; ωT]

## horizon
T = 21

## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

## initial conditions, controls, disturbances
ū = [0.01 * randn(m) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
storage = generate_storage(env.mechanism, x̄)
visualize(env.mechanism, storage; vis=vis)

# Objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-4 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(100.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

## problem 
prob = IterativeLQR.problem_data(model, obj)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# Solve
IterativeLQR.ilqr_solve!(prob,
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    verbose=true)

x_sol, u_sol = IterativeLQR.get_trajectory(prob)
storage = generate_storage(env.mechanism, x_sol)
visualize(env.mechanism, storage, vis=vis)
