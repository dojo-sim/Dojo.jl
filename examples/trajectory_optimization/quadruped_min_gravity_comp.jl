using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
include(joinpath(@__DIR__, "../../env/quadruped/methods/template.jl"))

gravity = -9.81
dt = 0.05
cf = 0.8 
damper = 50.0 
spring = 1.0
env = make("quadruped", 
    mode=:min, 
    dt=dt,
    g=gravity,
    cf=cf, 
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
d = 0 

# ## reference trajectory
initialize!(env.mechanism, :quadruped)
xref = quadruped_trajectory(env.mechanism, r=0.05, z=0.29; Δx=-0.04, Δfront=0.10, N=10, Ncycles=1)
zref = [min2max(env.mechanism, x) for x in xref]
visualize(env, xref)

# ## horizon 
T = 5

# ## model 
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
ū = [0.1 * randn(m) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū, w)
visualize(env, x̄)

# ## objective
qt = 1000.0 * ones(n)
ots = [(x, u, w) -> transpose(x - xref[1]) * Diagonal(qt) * (x - xref[1]) + transpose(u) * Diagonal(1.0e-3 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[1]) * Diagonal(qt) * (x - xref[1])

cts = IterativeLQR.Cost.(ots, n, m, d)
cT = IterativeLQR.Cost(oT, n, 0, 0)
obj = [cts..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - xref[1]
    return Δ
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
    verbose = true,
	linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0)

# ## solution
x_sol, u_sol = get_trajectory(prob)
visualize(env, x_sol)
