using IterativeLQR
using Dojo 

# ## system
dt = 0.05
gravity = -9.81
env = make("raiberthopper",
    mode=:min,
    dt=dt,
    g=gravity);

# ## visualizer
open(env.vis)

# ## dimensions
n = env.nx
m = env.nu
d = 0

# ## states
z1 = max2min(env.mechanism, raiberthopper_nominal_max())
zM = max2min(env.mechanism, raiberthopper_offset_max(0.5, 0.5, 1.0))
zT = max2min(env.mechanism, raiberthopper_offset_max(0.5, 0.5, 0.0))

# ## nominal control
u_control = [0.0; 0.0; env.mechanism.g * env.mechanism.timestep]

# ## horizon
T = 21

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w),
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
ū = [[0.0; 0.0; env.mechanism.g * env.mechanism.timestep + 0.0 * randn(1)[1]] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
visualize(env, x̄)

# ## objective
ot1 = (x, u, w) -> transpose(x - zM) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zM) + transpose(u) * Diagonal(0.1 * [0.1; 0.1; 0.1]) * u
ot2 = (x, u, w) -> transpose(x - zT) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zT) + transpose(u) * Diagonal(0.1 * [0.1; 0.1; 0.1]) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal([0.1; 0.1; 1.0; 0.001 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 1.0; 0.001]) * (x - zT)

ct1 = Cost(ot1, n, m, d)
ct2 = Cost(ot2, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - zT
    return  [Δ[collect(1:6)]; Δ[collect(12 .+ (1:2))]]
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
    con_tol=0.005,
    max_iter=100,
    max_al_iter=10,
    ρ_init=1.0,
    ρ_scale=10.0,
    verbose=true)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
visualize(env, x_sol)
