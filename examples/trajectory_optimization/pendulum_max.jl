using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
timestep = 0.05
gravity=-9.81
max_torque = 1.0e6
max_speed = 1.0e6
env = get_environment(:pendulum, 
    representation=:maximal, 
    max_speed=max_speed, 
    max_torque=max_torque,
    damper=0.0,
    timestep=timestep,
    gravity=gravity)

# ## visualizer 
open(env.vis)

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z1 = Dojo.pendulum_nominal_max()
zT = Dojo.pendulum_goal_max() 

# ## horizon
T = 26

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w, attitude_decompress=true),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w, attitude_decompress=true),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
ū = [10.0 * randn(m) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū)
visualize(env, x̄) 

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0 * ones(n)) * (x - zT)

ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints 
function ctrl_lmt(x, u, w) 
    [
    #  -max_torque - u[1]; 
    #  u[1] - max_torque;
    ]
end 

function goal(x, u, w)
    x - zT
end

# cont = IterativeLQR.Constraint(ctrl_lmt, n, m, 
#     idx_ineq=collect(1:2))
cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## solver
s = IterativeLQR.solver(model, obj, cons, 
    opts=IterativeLQR.Options(
        verbose=true,
        linesearch=:armijo,
        α_min=1.0e-5,
        obj_tol=1.0e-3,
        grad_tol=1.0e-3,
        max_iter=100,
        max_al_iter=10,
        ρ_init=1.0,
        ρ_scale=10.0))
IterativeLQR.initialize_controls!(s, ū)
IterativeLQR.initialize_states!(s, x̄)

# ## solve
IterativeLQR.solve!(s)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(s)

visualize(env, x_sol)


