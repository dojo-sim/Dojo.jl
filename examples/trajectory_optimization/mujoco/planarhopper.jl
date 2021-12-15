using Pkg; Pkg.activate(@__DIR__)
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCoViz 
using FiniteDiff

using IterativeLQR
using LinearAlgebra
using Random

# ## load MuJoCo model
path = joinpath(@__DIR__, "../../../env/raiberthopper/deps/planarhopper.xml")

include("mujoco_model.jl")
planarhopper_mujoco = MuJoCoModel(path)
sim = MJSim(planarhopper_mujoco.m, planarhopper_mujoco.d)

LyceumMuJoCoViz.visualize(sim) #, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

# ## horizon 
T = 11

# ## planarhopper_mujoco 
nx = planarhopper_mujoco.nx
nu = planarhopper_mujoco.nu 

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f!(y, planarhopper_mujoco, x, u), 
    (dx, x, u, w) -> fx!(dx, planarhopper_mujoco, x, u), 
    (du, x, u, w) -> fu!(du, planarhopper_mujoco, x, u), 
    nx, nx, nu) 

model = [dyn for t = 1:T-1] 

# ## initial conditions
x1 = []
xT = []

# ## objective
function objt(x, u, w)

	return J
end

function objT(x, u, w)
	return J
end

ct = IterativeLQR.Cost(objt, planarhopper_mujoco.nx, planarhopper_mujoco.nu, 0)
cT = IterativeLQR.Cost(objT, planarhopper_mujoco.nx, 0, 0)
obj = [[ct for t = 1:T-1]..., cT]

# ## constraints
goal(x, u, w) = x - xT

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, nx, 0)
cons = [[cont for t = 1:T-1]..., conT] 

# # rollout
Random.seed!(1)
ū = [1.0e-3 * randn(planarhopper_mujoco.nu) for t = 1:T-1]
w = [zeros(0) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū)

# ## problem
prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
IterativeLQR.reset!(prob.s_data)
IterativeLQR.solve!(prob, 
	linesearch = :armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-5,
    grad_tol=1.0e-5,
    max_iter=50,
    max_al_iter=10,
    con_tol=0.001,
    ρ_init=1.0, 
    ρ_scale=10.0,
	verbose=true)

@show prob.s_data.iter[1]
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)
@show prob.s_data.obj[1] # augmented Lagrangian cost


# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)

# ## MuJoCo visualizer
states = Array(undef, statespace(sim), T-1)
for t = 1:T-1
    sim.d.qpos .= x_sol[t][planarhopper_mujoco.idx_pos]
    sim.d.qvel .= x_sol[t][planarhopper_mujoco.idx_vel]
    sim.d.ctrl .= u_sol[t]
    states[:, t] .= getstate(sim)
end

visualize(sim, trajectories=[states]) # ctrl + LEFT (access trajectory mode)

