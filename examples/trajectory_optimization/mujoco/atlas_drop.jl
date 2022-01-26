using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5.xml")

include("mujoco_model.jl")
atlas = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(atlas.m, atlas.d)

# atlas.m.body_mass[1] += 1.0
# @show atlas.m.body_inertia[:, 1] = [0.01; 0.001; 0.001]
# length(atlas.d.ctrl)

T = 10
states = Array(undef, statespace(sim), T)
# sim.d.qpos .= [0.0; 1.17; 0.1; 0.1; -0.4; 0.4]
# sim.d.qvel .= zeros(6)
for t = 1:T
    # if t > 500 
    #     sim.d.qpos .= [1.0; 1.21; 0.0; 0.0; 0.0; 0.0] 
    # end
    # sim.d.qvel .= [0.0; 0.0; 0.0; 0.0]
    # sim.d.ctrl .= [0.0; 0.0; 0.1 * randn(1)[1]]#[5.0; 10.0; 0.0]
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
# sim.d.qpos
visualize(sim)#, trajectories=[states])

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
x1 = [0.0; 1.17; 0.1; 0.1; -0.4; 0.4; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xM = [0.5; 2.0; 0.1; 0.1; -0.4; 0.4; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
xT = [1.0; 1.17; 0.1; 0.1; -0.4; 0.4; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# # rollout
ū = [1.0 * randn(nu) for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective
obj1 = (x, u, w) -> (1 / T) * (transpose(x - xM) * Diagonal([1.0 * ones(6); 1.0 * ones(6)]) * (x - xM) + transpose(u) * Diagonal(1.0e-3 * [1.0; 1.0; 1.0]) * u)
obj2 = (x, u, w) -> (1 / T) * (transpose(x - xT) * Diagonal([1.0 * ones(6); 1.0 * ones(6)]) * (x - xT) + transpose(u) * Diagonal(1.0e-3 * [1.0; 1.0; 1.0]) * u)
objT = (x, u, w) -> transpose(x - xT) * Diagonal([1.0 * ones(6); 1.0 * ones(6)]) * (x - xT)

ct1 = IterativeLQR.Cost(obj1, nx, nu, 0)
ct2 = IterativeLQR.Cost(obj2, nx, nu, 0)
cT = IterativeLQR.Cost(objT, nx, 0, 0)
obj = [[ct1 for t = 1:500]..., [ct2 for t = 1:500]..., cT]

# ## constraints
goal(x, u, w) = (x - xT)

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
    max_iter=100,
    max_al_iter=3)

# ## solution
x_sol, u_sol = get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## 
using Plots
# plot(hcat(x_sol...)')
plot(hcat(u_sol..., u_sol[end])', linetype=:steppost)

# ## visualize
x_vis = [[x̄[1] for t = 1:10]..., x̄..., [x̄[end] for t = 1:10]...]
# x_vis = [[x_sol[1] for t = 1:100]..., x_sol..., [x_sol[end] for t = 1:100]...]
states = Array(undef, statespace(sim), length(x_vis))
for t = 1:length(x_vis)
    sim.d.qpos .= x_vis[t][1:6]
    sim.d.qvel .= x_vis[t][6 .+ (1:6)]
    # sim.d.ctrl .= 
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])

minimum([rank(fx) for fx in prob.m_data.model_deriv.fx])
minimum([rank(fu) for fu in prob.m_data.model_deriv.fu])