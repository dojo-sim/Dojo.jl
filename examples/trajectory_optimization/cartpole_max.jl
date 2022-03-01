
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra 

# ## system
gravity=-9.81
timestep = 0.1
env = get_environment(:cartpole, 
    representation=:maximal, 
    timestep=timestep,
    gravity=gravity);

## mujoco_inertia!(env.mechanism)

# ## visualizer 
open(env.vis) 

# ## dimensions
n = env.num_states
m = env.num_inputs

# ## states
z1 = cartpole_nominal_max()
zT = cartpole_goal_max()

# ## horizon
T = 26

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w), 
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m)
model = [dyn for t = 1:T-1]

# ## rollout
ū = [t < 5 ? 1.0 * rand(m) : (t < 10 ? -1.0 * rand(m) : zeros(m)) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū)
visualize(env, x̄)

# ## objective
ot = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT) + transpose(u) * Diagonal(1.0e-3 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(1.0e-1 * ones(n)) * (x - zT)

ct = IterativeLQR.Cost(ot, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
function goal(x, u, w)
    x - zT
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## solver 
solver = IterativeLQR.solver(model, obj, cons,
    opts=Options(
        max_al_iter=10,
        verbose=false))
IterativeLQR.initialize_controls!(solver, ū)
IterativeLQR.initialize_states!(solver, x̄)

# ## solve
@time IterativeLQR.solve!(solver)

# ## solution
z_sol, u_sol = IterativeLQR.get_trajectory(solver)
@show IterativeLQR.eval_obj(solver.m_data.obj.costs, solver.m_data.x, solver.m_data.u, solver.m_data.w)
@show solver.s_data.iter[1]
@show norm(goal(solver.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
visualize(env, [[z_sol[1] for t = 1:10]..., z_sol..., [z_sol[end] for t = 1:10]...])