# Open visualizer
vis = Visualizer()
open(vis)

using IterativeLQR

# System
dt = 0.05
gravity = -9.81
env = make("Hopper", 
    mode=:max, 
    dt=dt,
    g=gravity,
    vis=vis);

## state space
n = env.nx
m = env.nu
d = 0

function hopper_offset_max(x_shift, y_shift, z_shift)
    z = hopper_nominal_max()
    shift = [x_shift; y_shift; z_shift]
    z[1:3] += shift
    z[13 .+ (1:3)] += shift
    return z
end

z1 = hopper_nominal_max()
zM = hopper_offset_max(0.5, 0.5, 0.5)
zT = hopper_offset_max(0.5, 0.5, 0.0)

# nominal control
u_control = [0.0; 0.0; env.mechanism.g * env.mechanism.Δt]

# Time
T = 21
h = env.mechanism.Δt

# Model
function f(y, env::Environment, x, u, w)
	y .= step(env, x, u)[1]
end

function fx(fx, env::Environment, x, u, w)
	step(env, x, u, diff=true)
    fx .= env.fx
end

function fu(fu, env::Environment, x, u, w)
	step(env, x, u, diff=true)
	fu .= env.fu
end

dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# Initial conditions, controls, disturbances
ū = [[0.0; 0.0; env.mechanism.g * env.mechanism.Δt + 0.0 * randn(1)[1]] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
storage = generate_storage(env.mechanism, x̄)
visualize(env.mechanism, storage; vis=vis)

# Objective
ot1 = (x, u, w) -> transpose(x - zM) * Diagonal(dt * vcat([[0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)] for i=1:2]...)) * (x - zM) + transpose(u) * Diagonal(dt * [0.1; 0.1; 0.01]) * u
ot2 = (x, u, w) -> transpose(x - zT) * Diagonal(dt * vcat([[0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)] for i=1:2]...)) * (x - zT) + transpose(u) * Diagonal(dt * [0.1; 0.1; 0.01]) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(dt * vcat([[0.1 * ones(3); 0.001 * ones(3); 0.01 * ones(4); 0.01 * ones(3)] for i=1:2]...)) * (x - zT)

ct1 = Cost(ot1, n, m, d)
ct2 = Cost(ot2, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]

# Constraints
function goal(x, u, w)
    Δ = x - zT
    return [Δ[collect(1:6)]; Δ[collect(13 .+ (1:6))]]
end

cont = Constraint()
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# Solve
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

x_sol, u_sol = IterativeLQR.get_trajectory(prob)
storage = generate_storage(env.mechanism, x_sol)
visualize(env.mechanism, storage, vis=vis)
