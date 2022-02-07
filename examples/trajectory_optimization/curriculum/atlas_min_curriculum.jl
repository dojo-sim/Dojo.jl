using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
include(joinpath(module_dir(), "env/atlas/methods/template.jl"))

gravity = -9.81
dt = 0.05
cf = 0.8
damper = 50.0
spring = 0.0
ρ0 = 1e-2
env = atlas(
    mode=:min,
    dt=dt,
    gravity=gravity,
    cf=cf,
    damper=damper,
    spring=spring,
	model_type=:armless,
	infeasible_control=true,
	opts_step=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5),
    opts_grad=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5)
	)

env.mechanism.joints[2].constraints[2].damper = 100.0

# ## visualizer
open(env.vis)

# ## dimensions
n = env.nx
m = env.nu
d = 0

## simulate (test)
initialize!(env.mechanism, :atlas, tran=[1,0,0.0], rot=[0,0,0.], αhip=0.5, αknee=1.0)
function ctrl!(mech, k)
	u0 = -total_mass(env.mechanism) * env.mechanism.gravity* env.mechanism.timestep/1.1 * 0
	set_control!(mech, [u0; szeros(m-3)])
	return nothing
end
storage = simulate!(env.mechanism, 0.5, ctrl!, record=true, verbose=false,
	opts=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5))
visualize(env.mechanism, storage, vis=env.vis)

# ## reference trajectory
N = 1
initialize!(env.mechanism, :atlas, tran=[1,0,0.0], rot=[0,0,0.], αhip=0.5, αknee=1.0)
xref = atlas_trajectory(env.mechanism; timestep=dt, β=1.4, Δx=-0.03, r=0.08, z=0.89, N=10, Ncycles=N)
# # xref0 = deepcopy(xref[1])
# xref0[1:3] .+= [0.0, 0.0, 0.0] # body x
# xref0[4:6] .+= [0.0, 0.0, 0.0] # body ϕ
# # xref0[7:9] .+= [0.0, 0.0, 0.0] # body v
# # xref0[10:12] .+= [0.0, 0.0, 0.0] # body ϕ
# xref0[13:15] .+= [0.0, 0.0, 0.0] # rhip q
# # xref0[16:18] .+= [1.0, 0.0, 0.0] # rhip ϕ
# xref0[19:19] .+= [0.0] # rknee q
# # xref0[20:20] .+= [1.0] # rknee ϕ
# xref0[21:22] .+= [0.0, 0.0] # rankle q
# # xref0[23:24] .+= [0.0, 1.0] # rankle ϕ
# xref0[25:27] .+= [0.0, 0.0, 0.0] # lhip q
# # xref0[28:30] .+= [1.0, 1.0, 1.0] # lhip ϕ
# xref0[31:31] .+= [0.0] # lknee q
# # xref0[32:32] .+= [1.0] # lknee ϕ
# xref0[33:34] .+= [0.0, 0.0] # lankle q
# # xref0[35:36] .+= [0.0, 0.0] # lankle ϕ
# xref0[37:39] .+= [0.0, 0.0, 0.0] # pelvis q
# # xref0[40:42] .+= [1.0, 1.0, 1.0] # pelvis ϕ
# xref = [fill(xref[1], 3); fill(xref0, 10)]
visualize(env, xref)

# ## horizon
T = N * (21 - 1) + 1

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w),
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
u0 = -total_mass(env.mechanism) * env.mechanism.gravity* env.mechanism.timestep/1.5
ū = [[u0; zeros(m-3)] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū, w)
visualize(env, x̄)


# ## objective
qt = [
	[0.2, 0.2, 0.2, 0.05, 0.2, 0.2]; #x q body
	0.02ones(6); # v ϕ body
	0.8ones(3); 0.01ones(3); # rhip
	0.8ones(1); 0.01ones(1); # rknee
	0.8ones(2); 0.01ones(2); # rankle
	0.8ones(3); 0.01ones(3); # lhip
	0.8ones(1); 0.01ones(1); # lknee
	0.8ones(2); 0.01ones(2); # lankle
	0.5ones(3); 0.02ones(3); # pelvis
	]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(dt * qt) * (x - xref[t]) +
	transpose(u) * Diagonal(dt * [0.0001*ones(6); 0.02*ones(m-6)]) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(dt * qt) * (x - xref[end])

# ## constraints
cts = IterativeLQR.Cost.(ots, n, m, d)
cT = IterativeLQR.Cost(oT, n, 0, 0)
obj = [cts..., cT]


# ## constraints
function goal(x, u, w)
    Δ = x - xref[end]
    return Δ[[
		1,2,3,4,5,6, # body
		13,14,15, # rhip
		19, # rknee
		21,22, # rankle
		25,26,27, # lhip
		31, # lknee
		33,34, # lankle
		37,38,39 #pelvis
		]]
end

function ctrl_lmt(x, u, w)
	return u[collect(1:6)]
end

cont = IterativeLQR.Constraint(ctrl_lmt, n, m)
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## Initialization
u_prev = deepcopy(ū)
x_prev = deepcopy(xref) # TODO deepcopy(x̄)

# ## problem
prob = IterativeLQR.problem_data(model, obj, cons)

for ρ in [1e-2]#, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]
	println("ρ: ", scn(ρ), "   *************************************")
	IterativeLQR.initialize_controls!(prob, u_prev)
	IterativeLQR.initialize_states!(prob, x_prev)

	env.opts_step.rtol = ρ
	env.opts_step.btol = ρ
	env.opts_grad.rtol = ρ
	env.opts_grad.btol = ρ

	# ## solve
	@time IterativeLQR.solve!(prob,
	    verbose = true,
		linesearch=:armijo,
	    α_min=1.0e-5,
	    obj_tol=1.0e-3,
	    grad_tol=1.0e-3,
	    max_iter=100,
	    max_al_iter=5,
	    ρ_init=1.0,
	    ρ_scale=10.0)

	# Update initialization trajectory
	x_prev, u_prev = IterativeLQR.get_trajectory(prob)
end

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)
@show norm(vcat([ctrl_lmt(prob.m_data.x[t], prob.m_data.u[t], zeros(0)) for t=1:T-1]...), Inf)


# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(env.vis, cam_pos=[0,-3,2], zoom=3)


set_floor!(env.vis, x=3.0, y=1.0, z=0.01, alt=-0.015, color=RGBA(0,0.0,0.0,1.0))
set_camera!(env.vis, cam_pos=[0,-15,0], zoom=30)

convert_frames_to_video_and_gif("quadruped_clean_gait_side")
