using IterativeLQR

# ## system
include(joinpath(@__DIR__, "../../env/halfcheetah/methods/template.jl"))

dt = 0.05
gravity=-9.81
env = get_environment("halfcheetah", 
    mode=:minimal, 
    g=gravity,
    timestep=dt)

# ## visualizer 
open(env.vis)

# ## dimensions
n = env.num_states
m = env.num_inputs 
d = 0

# ## states
z1 = max2min(env.mechanism, halfcheetahState(x=0.00, z=0.00, θ=0.0))
zM = max2min(env.mechanism, halfcheetahState(x=0.25, z=0.40, θ=0.0))
zT = max2min(env.mechanism, halfcheetahState(x=0.50, z=0.00, θ=0.0))

# ## gravity compensation control
mech = getmechanism(:halfcheetah, timestep=dt, g=gravity, damper=100.0, spring=1000.0)
initialize!(mech, :halfcheetah, x=0.0, z=0.0, θ=0.0)
storage = simulate!(mech, 2.0, record=true, verbose=false)
visualize(mech, storage, vis=env.vis)
ugc = gravity_compensation(mech)
u_control = ugc[3 .+ (1:6)]

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
ū = [u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, z1, ū, w)
visualize(env, x̄) 

# ## objective
qt1 = [0.1; 0.1; 1.0; 0.01 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 0.001 * (6)]
qt2 = [0.1; 0.1; 1.0; 0.01 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 0.001 * (6)]
qt  = [0.1; 0.1; 1.0; 0.01 * ones(3); 0.001 * ones(3); 0.01 * ones(3); 0.001 * (6)]

ot1 = (x, u, w) -> transpose(x - zM) * Diagonal(dt * qt) * (x - zM) + transpose(u) * Diagonal(dt * 0.01 * ones(m)) * u
ot2 = (x, u, w) -> transpose(x - zT) * Diagonal(dt * qt) * (x - zT) + transpose(u) * Diagonal(dt * 0.01 * ones(m)) * u
oT = (x, u, w) -> transpose(x - zT) * Diagonal(dt * qt) * (x - zT)

ct1 = Cost(ot1, n, m, d)
ct2 = Cost(ot2, n, m, d)
cT = Cost(oT, n, 0, 0)
obj = [[ct1 for t = 1:10]..., [ct2 for t = 1:10]..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - zT
    return Δ[collect(1:6)]
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
    verbose=true,
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
