using IterativeLQR

# ## system
include(joinpath(@__DIR__, "../../env/quadruped/methods/template.jl"))

gravity = -9.81
dt = 0.05
cf = 0.8 
damper = 5.0 
spring = 0.0
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
initialize!(env.mechanism, :quadruped)
storage = simulate!(env.mechanism, 0.5, record=true, verbose=false)
visualize(env.mechanism, storage, vis=env.vis)

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
T = 21 

# ## model 
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w), 
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
ū = [zeros(m) for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = rollout(model, x1, ū, w)
visualize(env, x̄)

# Objective
qt = [0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill([0.2, 0.001], 12)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(dt * qt) * (x - xref[t]) + transpose(u) * Diagonal(dt * 0.01 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(dt * qt) * (x - xref[end])

cts = Cost.(ots, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [cts..., cT]

# Constraints
function goal(x, u, w)
    Δ = x - xref[end]
    return Δ[collect(1:3)]
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

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

x_sol, u_sol = get_trajectory(prob)
visualize(env, x_sol)
