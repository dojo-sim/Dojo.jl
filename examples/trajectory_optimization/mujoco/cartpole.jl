using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/cartpole/deps/cartpole.xml")

include("mujoco_model.jl")
cartpole = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(cartpole.m, cartpole.d)

@show cartpole.m.body_mass 
@show cartpole.m.body_inertia[:, :]

LyceumMuJoCoViz.visualize(sim) #, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

# ## horizon 
T = 26

# ## cartpole 
nx = cartpole.nx
nu = cartpole.nu 

# ## model
dyn = Dynamics(
    (y, x, u, w) -> f!(y, cartpole, x, u), 
    (dx, x, u, w) -> fx!(dx, cartpole, x, u), 
    (du, x, u, w) -> fu!(du, cartpole, x, u), 
    nx, nx, nu) 

model = [dyn for t = 1:T-1] 

# ## initial conditions
x1 = [0.0; π; 0.0; 0.0]
xM = zeros(cartpole.nx)
xT = zeros(cartpole.nx)
x_ref = zeros(cartpole.nx)

# x1 = [0.0; 0.0; 0.0; 0.0]
# xM = [0.0; 0.2 * π; 0.0; 0.0]
# xT = [0.0; 0.2 * π; 0.0; 0.0]
# x_ref = [0.0; 0.2 * π; 0.0; 0.0]

# # rollout
ū = [t < 5 ? 1.0 * rand(nu) : (t < 10 ? -1.0 * rand(nu) : zeros(nu)) for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = rollout(model, x1, ū)

# ## objective

function objt(x, u, w)
	J = 0.0 
	J += transpose(x - x_ref) * Diagonal(1.0e-1 * ones(nx)) * (x - x_ref) 
	J += transpose(u) * Diagonal(1.0e-3 * ones(nu)) * u
	return J
end

function objT(x, u, w)
	J = 0.0 
	J += transpose(x - x_ref) * Diagonal(1.0e-1 * ones(nx)) * (x - x_ref) 
	return J
end

ct = IterativeLQR.Cost(objt, nx, nu, 0)
cT = IterativeLQR.Cost(objT, nx, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
goal(x, u, w) = x - xT

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
    max_al_iter=10)

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
    sim.d.qpos .= x_vis[t][1:2]
    sim.d.qvel .= x_vis[t][3:4]
    sim.d.ctrl .= [0.0]
    # LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
end
visualize(sim, trajectories=[states])