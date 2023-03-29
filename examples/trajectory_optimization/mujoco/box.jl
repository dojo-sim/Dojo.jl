using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/box/deps/block.xml")

include("mujoco_model.jl")
box = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(box.m, box.d)

@show box.m.body_mass 
@show box.m.body_inertia[:, :]

T = 100
states = Array(undef, statespace(sim), T)
sim.d.qpos[:] = [0.0; 0.0; -0.24 + 1.0; 0.0; 0.0; 0.0] 
sim.d.qvel[:] = zeros(6)
sim.d.ctrl
for t = 1:T
    sim.d.ctrl[:] = [1.0; 0.0; 0.0]
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])

# ## horizon 
T = 101

# ## box 
nx = box.nx
nu = box.nu 

# ## model
dyn = Dynamics(
    (y, x, u, w) -> f!(y, box, x, u), 
    (dx, x, u, w) -> fx!(dx, box, x, u), 
    (du, x, u, w) -> fu!(du, box, x, u), 
    nx, nx, nu) 

model = [dyn for t = 1:T-1] 

# ## initial conditions
x1 = [0.0; 0.0; -0.24; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
# xT = [1.0; 0.0; -0.24; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # right goal
xT = [0.0; 0.0; -0.24 + 1.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0] # right goal

# x1 = [0.0; 0.0; 0.0; 0.0]
# xM = [0.0; 0.2 * π; 0.0; 0.0]
# xT = [0.0; 0.2 * π; 0.0; 0.0]
# x_ref = [0.0; 0.2 * π; 0.0; 0.0]

# # rollout
ū = [zeros(3) for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective
ot = (x, u, w) -> transpose(x - xT) * Diagonal(1.0 * ones(nx)) * (x - xT) + transpose(u) * Diagonal(1.0e-2 * ones(nu)) * u
oT = (x, u, w) -> transpose(x - xT) * Diagonal(1.0 * ones(nx)) * (x - xT)

ct = Cost(ot, nx, nu)
cT = Cost(oT, nx, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
goal(x, u, w) = x - xT

cont = Constraint()
conT = Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# ## problem
prob = solver(model, obj, cons, 
    opts=Options(
        linesearch=:armijo,
        α_min=1.0e-5,
        obj_tol=1.0e-3,
        grad_tol=1.0e-3,
        con_tol=0.005,
        max_iter=100,
        max_al_iter=10,
        ρ_init=1.0,
        ρ_scale=10.0,
        verbose=true))
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
@time IterativeLQR.solve!(prob)

# ## solution
x_sol, u_sol = get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## 
# using Plots
# plot(hcat(x_sol...)')
# plot(hcat(u_sol..., u_sol[end])', linetype=:steppost)

# ## visualize
x_vis = [[x_sol[1] for t = 1:10]..., x_sol..., [x_sol[end] for t = 1:10]...]
states = Array(undef, statespace(sim), length(x_vis))
for t = 1:length(x_vis)
    sim.d.qpos .= x_vis[t][1:6]
    sim.d.qvel .= x_vis[t][6 .+ (1:6)]
    sim.d.ctrl .= [0.0]
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])