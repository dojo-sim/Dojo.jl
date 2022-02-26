using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
include(joinpath(module_dir(), "environments/atlas/methods/template.jl"))

gravity=-9.81
dt = 0.05
friction_coefficient = 0.8
damper = 20.0
spring = 0.0
ρ0 = 3e-1
env = atlas(
    mode=:minimal,
    timestep=dt,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
	model_type=:armless,
	infeasible_control=true,
	opts_step=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5),
    opts_grad=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5)
	)

env.mechanism.joints[2].rotational.damper = 75.0
env.mechanism.joints[7].rotational.damper = 75.0
env.mechanism.joints[8].rotational.damper = 75.0

# ## visualizer
open(env.vis)

# ## dimensions
n = env.num_states
m = env.num_inputs
d = 0

## simulate (test)
initialize!(env.mechanism, :atlas, model_type=:armless, tran=[1,0,0.0], rot=[0,0,0.], αhip=0.5, αknee=1.0)
function ctrl!(mech, k)
	u0 = -total_mass(env.mechanism) * env.mechanism.gravity* env.mechanism.timestep/1.1 * 0
	set_input!(mech, [u0; szeros(m-3)])
	return nothing
end
storage = simulate!(env.mechanism, 0.5, ctrl!, record=true, verbose=false,
	opts=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5))
visualize(env.mechanism, storage, vis=env.vis, show_contact=false)

mech.contacts[1].model
# ## reference trajectory
N = 3
initialize!(env.mechanism, :atlas, model_type=:armless,
	tran=[0,0,0.0], rot=[0,0,0.], αhip=0.5, αknee=1.0)
xref = atlas_trajectory(env.mechanism; timestep=dt, β=1.4,
	αtorso=0.07, Δx=-0.03, r=0.08, z=1.12, N=10, Ncycles=N)
# x = get_minimal_state(env.mechanism)
# xref0 = deepcopy(xref[1])
# # xref0[1:3] .+= [1.0, 1.0, 1.0] # floating x
# # xref0[4:6] .+= [1.0, 0.0, 1.0] # floating ϕ
# # xref0[7:9] .+= [0.0, 0.0, 0.0] # floating v
# # xref0[10:12] .+= [0.0, 0.0, 0.0] # floating ϕ
# xref0[13:15] .+= [0.0, 0.0, 0.0] # back q
# # xref0[16:18] .+= [1.0, 0.0, 0.0] # back ϕ
# xref0[19:21] .+= [0.0, 1.0, 0.0] # rhip q
# # xref0[22:24] .+= [1.0, 0.0, 0.0] # rhip ϕ
# xref0[25:25] .+= [0.0] # rknee q
# # xref0[26:26] .+= [1.0] # rknee ϕ
# xref0[27:28] .+= [0.0, 0.0] # rankle q
# # xref0[29:30] .+= [0.0, 1.0] # rankle ϕ
# xref0[31:33] .+= [0.0, 0.0, 0.0] # lhip q
# # xref0[34:36] .+= [1.0, 1.0, 1.0] # lhip ϕ
# xref0[37:37] .+= [0.0] # lknee q
# # xref0[38:38] .+= [1.0] # lknee ϕ
# xref0[39:40] .+= [0.0, 0.0] # lankle q
# # xref0[41:42] .+= [0.0, 0.0] # lankle ϕ
# xref0 = deepcopy(x)
# xref0[13] += 1.0
# xref = [fill(xref[1], 3); fill(xref0, 10)]
# visualize(env, fill(xref0, 10))
visualize(env, xref)

# ## horizon
T = N * (21 - 1) + 1
# T = 51

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w),
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
u0 = -total_mass(env.mechanism) * env.mechanism.gravity* env.mechanism.timestep/1.1
ū = [[u0; 0; -0.4; 0; zeros(m-6)] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū, w)
visualize(env, x̄)


# ## objective
qt = [
	[0.2, 0.05, 0.2, 0.05, 0.2, 0.2]; #x q floating
	0.02ones(6); # v ϕ floating
	0.5ones(3); 0.02ones(3); # back
	4.8ones(3); 0.01ones(3); # rhip
	4.8ones(1); 0.01ones(1); # rknee
	4.8ones(2); 0.01ones(2); # rankle
	4.8ones(3); 0.01ones(3); # lhip
	4.8ones(1); 0.01ones(1); # lknee
	4.8ones(2); 0.01ones(2); # lankle
	]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(dt * qt) * (x - xref[t]) +
	transpose(u) * Diagonal(dt * [0.01*ones(6); 0.02*ones(m-6)]) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(dt * qt) * (x - xref[end])

# ## constraints
cts = IterativeLQR.Cost.(ots, n, m, d)
cT = IterativeLQR.Cost(oT, n, 0, 0)
obj = [cts..., cT]


# ## constraints
function goal(x, u, w)
    Δ = x - xref[end]
    return Δ[[
		1,2,3,4,5,6, # floating
		13,14,15, # back
		19,20,21, #pelvis
		25, # rknee
		27,28, # rankle
		31,32,33, # lhip
		37, # lknee
		39,40, # lankle
		]]
end

function ctrl_lmt(x, u, w)
	return 1e-1*u[collect(1:6)]
end

cont = IterativeLQR.Constraint(ctrl_lmt, n, m)
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## Initialization
u_prev = deepcopy(ū)
# x_prev = deepcopy(x̄)
x_prev = deepcopy(xref) # TODO deepcopy(x̄)

# ## problem
prob = IterativeLQR.problem_data(model, obj, cons)

for ρ in [1e-3, 3e-4, 1e-4]#, 3e-5]#, 1e-5, 3e-6, 1e-6]
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
		max_al_iter=4,
	    ρ_init=1.0,
	    ρ_scale=10.0)

	# Update initialization trajectory
	x_prev, u_prev = IterativeLQR.get_trajectory(prob)
end
scatter(abs.(goal(prob.m_data.x[T], zeros(0), zeros(0))))
# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)
@show norm(vcat([ctrl_lmt(prob.m_data.x[t], prob.m_data.u[t], zeros(0)) for t=1:T-1]...), Inf)

# jldsave(joinpath(@__DIR__, "atlas_traj_6steps.jld2"), x_sol=x_sol, u_sol=u_sol)

plot([x[5] for x in x_sol])


# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(env.vis, cam_pos=[0,-3,2], zoom=3)
open(env.vis)

set_floor!(env.vis, x=6.0, y=6.0, z=0.01, alt=0.0, color=RGBA(0.5,0.5,0.5,1.0))
set_camera!(env.vis, cam_pos=[5,-6,6], zoom=3)
set_light!(env.vis)

convert_frames_to_video_and_gif("atlas_6_steps_side")

render_static(env.vis)
# open(joinpath(@__DIR__, "atlas_6_steps.html"), "w") do file
#     write(file, static_html(env.vis))
# end
