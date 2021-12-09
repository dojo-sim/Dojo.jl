# Open visualizer
vis = Visualizer()
open(vis)

using IterativeLQR

# system
dt = 0.05
gravity = -9.81
max_torque = 20.0
env = make("pendulum", 
    mode=:min, 
    dt=dt,
    g=gravity,
    max_torque=max_torque,
    vis=vis);

## state space
n = env.nx
m = env.nu
d = 0

## initial state
z1 = [0.0; 0.0]

## goal state 
zT = [π; 0.0]

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
ū = [0.1 * randn(m) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
visualize(env, x̄) 

# Objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-3 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT)

ct = Cost(ot, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# Constraints
function goal(x, u, w)
    x - zT
end

cont = Constraint()
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

## problem 
prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# Solve
IterativeLQR.solve!(prob,
    verbose = true,
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0)

x_sol, u_sol = IterativeLQR.get_trajectory(prob)
visualize(env, x_sol)