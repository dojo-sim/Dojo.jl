using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
include(joinpath(@__DIR__, "../../env/quadruped/methods/template.jl"))

gravity = -9.81
dt = 0.05
friction_coefficient = 0.8
damper = 5.0
spring = 0.0
env = quadruped(
    mode=:min,
    dt=dt,
    body_contact=false,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring)

# ## visualizer
open(env.vis)

# ## simulate (test)
# initialize!(env.mechanism, :quadruped)
# storage = simulate!(env.mechanism, 0.5, record=true, verbose=false)
# visualize(env.mechanism, storage, vis=env.vis)

# ## dimensions
n = env.nx
m = env.nu

# ## reference trajectory
N = 2
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism, r=0.05, z=0.29; Δx=-0.04, Δfront=0.10, N=10, Ncycles=N)
zref = [minimal_to_maximal(env.mechanism, x) for x in xref]
visualize(env, xref)

## gravity compensation TODO: solve optimization problem instead
mech = get_mechanism(:quadruped, timestep=dt, gravity=gravity, friction_coefficient=0.8, damper=5.0, spring=0.0)
initialize!(mech, :quadruped)
storage = simulate!(mech, 1.0, record=true, verbose=false)
visualize(mech, storage, vis=env.vis)
ugc = gravity_compensation(mech)
u_control = ugc[6 .+ (1:12)]

# ## horizon
T = N * (21 - 1) + 1

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w),
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
ū = [u_control for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū)
visualize(env, x̄)

# ## objective
qt = [0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(dt * qt) * (x - xref[t]) + transpose(u) * Diagonal(dt * 0.01 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(dt * qt) * (x - xref[end])

cts = IterativeLQR.Cost.(ots, n, m)
cT = IterativeLQR.Cost(oT, n, 0)
obj = [cts..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - xref[end]
    return Δ[collect(1:3)]
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem
prob = IterativeLQR.solver(model, obj, cons, 
    opts=Options(verbose = false,
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
@time IterativeLQR.solve!(prob)

vis = Visualizer()
open(env.vis)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(env.vis, cam_pos=[0,-3,2], zoom=3)


x_shift = deepcopy(x_sol)
for x in x_shift
    x[3] += 0.01
end
z = [minimal_to_maximal(env.mechanism, x) for x in x_shift]

t = 1 #10, 20, 30, 41
set_robot(env.vis, env.mechanism, z[41])

# ## visualize
x_view = [[x_shift[1] for t = 1:15]..., x_shift..., [x_shift[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(vis, cam_pos=[0,-50,0], zoom=30)
set_floor!(env.vis, z=0.01)
