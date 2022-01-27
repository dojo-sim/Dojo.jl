using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
using JLD2
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

# @show hopper.m.body_mass[2] = 1.0 
# @show hopper.m.body_mass[3] = 0.1 
# @show hopper.m.body_inertia[:, 2] .= 0.004 
# @show hopper.m.body_inertia[:, 3] .= 0.0001

@show hopper.m.body_mass[:]
@show hopper.m.body_inertia[:, :]

# @show length(hopper.d.ctrl)

# LyceumMuJoCoViz.visualize(sim) #, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

T = 100
states = Array(undef, statespace(sim), T)
sim.d.qpos[:] = [0.0; 0.0; 0.55; 0.0; 0.0; 0.0; 0.5] 
sim.d.qvel[:] = zeros(7)
sim.d.ctrl
for t = 1:T
    sim.d.ctrl[:] = [0.0; 0.0; 1.0 * hopper.m.body_mass[2] * 9.81 / 200; 0.0; 0.0; 0.0; 0.0]
    # sim.d.ctrl .= 1.0 * randn(3)#[0.0; 0.0; 1.0 * hopper.m.body_mass[2] * 9.81]#[5.0; 10.0; 0.0]
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])
# 2.0 * hopper.m.body_mass[2] * 9.81

# ## horizon 
T = 101
Tm = convert(Int, floor((T - 1) / 2))

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
xM = [0.5; 0.5; 0.55 + 0.5; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xT = [0.5; 0.5; 0.55 + 0.0; 0.0; 0.0; 0.0; 0.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# # rollout
ū = [[0.0; 0.0; 1.0 * hopper.m.body_mass[2] * 9.81 / 200; 0.0; 0.0; 0.0; 0.0] for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective
obj1 = (x, u, w) -> (transpose(x - xM) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001]) * (x - xM) + transpose(u) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]) * u)
obj2 = (x, u, w) -> (transpose(x - xT) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001]) * (x - xT) + transpose(u) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]) * u)
objT = (x, u, w) -> (transpose(x - xT) * Diagonal([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001]) * (x - xT))

ct1 = IterativeLQR.Cost(obj1, nx, nu, 0)
ct2 = IterativeLQR.Cost(obj2, nx, nu, 0)
cT = IterativeLQR.Cost(objT, nx, 0, 0)
obj = [[ct1 for t = 1:Tm]..., [ct2 for t = 1:Tm]..., cT]

# ## constraints
god_control(x, u, w) = u[4:7]
goal(x, u, w) = 1.0 * (x - xT)[collect(1:7)]

jcont = Constraint(god_control, nx, nu)
conT = Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# ## problem
prob = problem_data(model, obj, cons)
initialize_controls!(prob, ū)
initialize_states!(prob, x̄)

# ## solve
@time solve!(prob, 
    linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    con_tol=0.005,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0, 
    verbose=true)

# ## solution
x_sol, u_sol = get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

path = joinpath(@__DIR__, "hopper_solution.jld2")
# @save path x_sol u_sol

# ## 
# using Plot
# plot(hcat(x_sol...)')
# plot(hcat(u_sol..., u_sol[end])', linetype=:steppost)

# ## visualize
x_vis = [[x̄[1] for t = 1:10]..., x̄..., [x̄[end] for t = 1:10]...]
x_vis = [[x_sol[1] for t = 1:100]..., x_sol..., [x_sol[end] for t = 1:100]...]
states = Array(undef, statespace(sim), length(x_vis))
for t = 1:length(x_vis)
    sim.d.qpos[:] .= x_vis[t][1:7]
    sim.d.qvel[:] .= x_vis[t][7 .+ (1:7)]
    # sim.d.ctrl .= 
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])

using Plots
plot(hcat(u_sol...)[1:3, :]')
plot(hcat(u_sol...)[4:7, :]')