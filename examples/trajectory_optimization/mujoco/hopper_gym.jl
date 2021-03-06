using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/hopper/deps/hopper.xml")

include("mujoco_model.jl")
hopper = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(hopper.m, hopper.d)

# @show hopper.m.body_mass 
# @show hopper.m.body_inertia[:, :]
# length(hopper.d.ctrl)

# LyceumMuJoCoViz.visualize(sim) #, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

T = 100
states = Array(undef, statespace(sim), T)
sim.d.qpos #.= 
sim.d.qvel #.= zeros(7)
sim.d.ctrl
for t = 1:T
    # sim.d.qpos .= [0.0; 0.0; 0.55; 0.5]
    # sim.d.qvel .= [0.0; 0.0; 0.0; 0.0]
    sim.d.ctrl .= 1.0 * randn(3)#[0.0; 0.0; 1.0 * hopper.m.body_mass[2] * 9.81]#[5.0; 10.0; 0.0]
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])
# 2.0 * hopper.m.body_mass[2] * 9.81

# ## horizon 
T = 201

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
x1 = [0.0; 1.21; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xM = [0.25; 2.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]#[0.5; 0.5; 1.05; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xT = [0.5; 1.21; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]#[0.5; 0.5; 0.55; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# # rollout
ū = [0.01 * rand(3) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective
obj1 = (x, u, w) -> transpose(x - xM) * Diagonal([10.0 * ones(6); 1.0e-2 * ones(6)]) * (x - xM) + transpose(u) * Diagonal([1.0; 1.0; 1.0]) * u
obj2 = (x, u, w) -> transpose(x - xT) * Diagonal([10.0 * ones(6); 1.0e-2 * ones(6)]) * (x - xT) + transpose(u) * Diagonal([1.0; 1.0; 1.0]) * u
objT = (x, u, w) -> transpose(x - xT) * Diagonal([ones(6); ones(6)]) * (x - xT)

ct1 = IterativeLQR.Cost(obj1, nx, nu)
ct2 = IterativeLQR.Cost(obj2, nx, nu)
cT = IterativeLQR.Cost(objT, nx, 0)
obj = [[ct1 for t = 1:100]..., [ct2 for t = 1:100]..., cT]

# ## constraints
goal(x, u, w) = (x - xT)[collect(1:6)]

cont = Constraint()
conT = Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# ## problem
prob = solver(model, obj, cons, 
    opts=Options(
        verbose=true,
        max_al_iter=5))
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
solve!(prob)

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
x_vis = [[x̄[1] for t = 1:10]..., x̄..., [x̄[end] for t = 1:10]...]
x_vis = [[x_sol[1] for t = 1:100]..., x_sol..., [x_sol[end] for t = 1:100]...]
states = Array(undef, statespace(sim), length(x_vis))
for t = 1:length(x_vis)
    sim.d.qpos .= x_vis[t][1:6]
    sim.d.qvel .= x_vis[t][6 .+ (1:6)]
    # sim.d.ctrl .= 
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])