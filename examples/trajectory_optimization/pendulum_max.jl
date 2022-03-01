using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
timestep = 0.05
gravity=-9.81
max_torque = 200.0
max_speed = 8.0
env = get_environment(:pendulum, 
    representation=:maximal, 
    max_speed=max_speed, 
    max_torque=max_torque,
    damper=1.0,
    timestep=timestep,
    gravity=gravity,
    vis=vis);

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z1 = pendulum_nominal_max()
zT = pendulum_goal_max() 

# ## horizon
T = 51

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
ū = [1.0 * randn(m) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū)
visualize(env, x̄) 

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-2 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT)

ct = Cost(ot, n, m)
cT = Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints 
function ctrl_lmt(x, u, w) 
    [
     -max_torque - u[1]; 
     u[1] - max_torque;
    ]
end 

function goal(x, u, w)
    x - zT
end

cont = Constraint(ctrl_lmt, n, m, idx_ineq=collect(1:2))
conT = Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem 
prob = IterativeLQR.solver(model, obj, cons, 
    opts=Options(verbose = true,
        linesearch=:armijo,
        α_min=1.0e-5,
        obj_tol=1.0e-3,
        grad_tol=1.0e-3,
        max_iter=100,
        max_al_iter=5,
        ρ_init=1.0,
        ρ_scale=10.0))
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
IterativeLQR.solve!(prob)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
visualize(env, x_sol)


