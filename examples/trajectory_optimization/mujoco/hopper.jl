using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/raiberthopper/deps/hopper.xml")

include("mujoco_model.jl")
hopper = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(hopper.m, hopper.d)

@show hopper.m.body_mass 
@show hopper.m.body_inertia[:, :]
length(hopper.d.ctrl)

# LyceumMuJoCoViz.visualize(sim) #, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

T = 1000
states = Array(undef, statespace(sim), T)
sim.d.qpos #.= [0.0; 0.0; 0.55; 0.0; 0.0; 0.0; 0.5]
sim.d.qvel .= zeros(7)
for t = 1:T
    # sim.d.qpos .= [0.0; 0.0; 0.55; 0.5]
    # sim.d.qvel .= [0.0; 0.0; 0.0; 0.0]
    # sim.d.ctrl .= [0.0; 0.0; 1.0 * hopper.m.body_mass[2] * 9.81]#[5.0; 10.0; 0.0]
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])
2.0 * hopper.m.body_mass[2] * 9.81

# ## horizon 
T = 1001

# ## hopper 
nx = hopper.nx
nu = hopper.nu 

# ## model
dyn = Dynamics(
    (y, x, u, w) -> f!(y, hopper, x, u), 
    (dx, x, u, w) -> fx!(dx, hopper, x, u), 
    (du, x, u, w) -> fu!(du, hopper, x, u), 
    nx, nx, nu) 

model = [dyn for t = 1:T-1] 

# ## initial conditions
x1 = [0.0; 0.0; 0.55; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xM = [0.1; 0.0; 0.55 + 0.25; 0.0; 0.0; 0.0; 0.25; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]#[0.5; 0.5; 1.05; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xT = [0.1; 0.0; 0.55 + 0.25; 0.0; 0.0; 0.0; 0.25; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]#[0.5; 0.5; 0.55; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# # rollout
ū = [[0.0; 0.0; hopper.m.body_mass[2] * 9.81] for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective
dt = 0.001
obj1 = (x, u, w) -> transpose(x - xM) * Diagonal(dt * [0.1; 0.1; 1.0; 0.01 * ones(3); 1.0; 0.01 * ones(3); 0.01 * ones(3); 0.001]) * (x - xM) + transpose(u) * Diagonal(dt * [0.01; 0.01; 0.01]) * u
obj2 = (x, u, w) -> transpose(x - xT) * Diagonal(dt * [0.1; 0.1; 1.0; 0.01 * ones(3); 1.0; 0.01 * ones(3); 0.01 * ones(3); 0.001]) * (x - xT) + transpose(u) * Diagonal(dt * [0.01; 0.01; 0.01]) * u
objT = (x, u, w) -> transpose(x - xT) * Diagonal(dt * [0.1; 0.1; 1.0; 0.01 * ones(3); 1.0; 0.01 * ones(3); 0.01 * ones(3); 0.001]) * (x - xT)

ct1 = IterativeLQR.Cost(obj1, nx, nu, 0)
ct2 = IterativeLQR.Cost(objt, nx, nu, 0)
cT = IterativeLQR.Cost(objT, nx, 0, 0)
obj = [[ct1 for t = 1:500]..., [ct2 for t = 1:500]..., cT]

# ## constraints
goal(x, u, w) = (x - xT)[collect(1:7)]

cont = Constraint()
conT = Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# ## problem
prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
solve!(prob, 
    verbose=true,
    max_al_iter=7)

# ## solution
x_sol, u_sol = get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## 
# using Plot
# plot(hcat(x_sol...)')
# plot(hcat(u_sol..., u_sol[end])', linetype=:steppost)

# ## visualize
# x_vis = [[x̄[1] for t = 1:10]..., x̄..., [x̄[end] for t = 1:10]...]
x_vis = [[x_sol[1] for t = 1:100]..., x_sol..., [x_sol[end] for t = 1:100]...]
states = Array(undef, statespace(sim), length(x_vis))
for t = 1:length(x_vis)
    sim.d.qpos .= x_vis[t][1:7]
    sim.d.qvel .= x_vis[t][7 .+ (1:7)]
    # sim.d.ctrl .= 
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])