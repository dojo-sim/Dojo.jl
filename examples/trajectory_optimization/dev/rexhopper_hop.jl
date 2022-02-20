using Dojo
using IterativeLQR
using LinearAlgebra

# ## system
# include(joinpath(module_dir(), "env/atlas/methods/template.jl"))
include(joinpath(module_dir(), "env/rexhopper/methods/env.jl"))
include(joinpath(module_dir(), "env/rexhopper/methods/initialize.jl"))


gravity = -9.81
dt = 0.02
friction_coefficient = 0.5
damper = 2.0
spring = 5.0
ρ0 = 1e-4
env = rexhopper(
    mode=:min,
    timestep=dt,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
	contact=true,
	contact_body=true,
	# model=:rexhopper_fixed,
	model=:rexhopper_no_wheel,
	infeasible_control=true,
	opts_step=SolverOptions(rtol=ρ0, btol=ρ0, undercut=5.0),
    opts_grad=SolverOptions(rtol=ρ0, btol=ρ0, undercut=5.0)
	)


# ## visualizer
open(env.vis)

# ## dimensions
n = env.nx
m = env.nu
d = 0

## simulate (test)
initialize!(env.mechanism, :rexhopper, x=[0,0,0], θ=[0,0,0.])
xinit = get_minimal_state(env.mechanism) + 0.3*[ones(6); zeros(16)]# 0.1ones(minimal_dimension(env.mechanism))
zinit = get_maximal_state(env.mechanism)
x2 = maximal_to_minimal(env.mechanism, zinit)
z2 = minimal_to_maximal(env.mechanism, xinit)
# build_robot(env.vis, env.mechanism)
set_robot(env.vis, env.mechanism, zinit)
set_robot(env.vis, env.mechanism, z2)


initialize!(env.mechanism, :rexhopper, x=[0,0,0.0], θ=[0,0,0.])
xup = get_minimal_state(env.mechanism)
function ctrl!(mech, k)
	u0 = -total_mass(env.mechanism) * env.mechanism.gravity* env.mechanism.timestep/1.1 * 0
	set_control!(mech, [u0; szeros(m-3)])
	return nothing
end
storage = simulate!(env.mechanism, 0.5, ctrl!, record=true, verbose=false,
	opts=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5))
visualize(env.mechanism, storage, vis=env.vis, show_contact=false)

# ## horizon
T = 25
# ## reference trajectory
xref = [deepcopy(xup) for i = 1:T]
visualize(env, xref)

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> f(y, env, x, u, w),
    (dx, x, u, w) -> fx(dx, env, x, u, w),
    (du, x, u, w) -> fu(du, env, x, u, w),
    n, n, m)

model = [dyn for t = 1:T-1]

# ## rollout
initialize!(env.mechanism, :rexhopper, x=[0,0,0.0], θ=[0,0,0.])
xinit = get_minimal_state(env.mechanism)
x1 = xinit
u0 = -total_mass(env.mechanism) * env.mechanism.gravity * env.mechanism.timestep/1.1
ū = [[u0; 0; 0.5env.mechanism.timestep ; 0; zeros(m-6)] for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū, w)
visualize(env, x̄)


# ## objective
qt = 0.01*[
	[0.2, 0.05, 0.2, 0.05, 0.2, 0.2]; #x q floating
	0.02ones(6); # v ϕ floating
	0.8ones(0); 0.01ones(0); # joint_rwz
	0.8ones(0); 0.01ones(0); # joint_rw1
	0.8ones(0); 0.01ones(0); # joint_rw0
	0.8ones(1); 0.01ones(1); # joint2
	0.8ones(1); 0.01ones(1); # joint3
	0.8ones(1); 0.01ones(1); # joint0
	0.8ones(1); 0.01ones(1); # joint1
	0.8ones(1); 0.01ones(1); # loop_joint
	]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(dt * qt) * (x - xref[t]) +
	transpose(u) * Diagonal(dt * [0.1*ones(6); 0.2*ones(m-6)]) * u for t = 1:T-1]
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
		]]
end

function ctrl_lmt(x, u, w)
	return 1e-1*u[collect(1:3)]
end

# cont = IterativeLQR.Constraint(ctrl_lmt, n, m)
cont = Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## Initialization
u_prev = deepcopy(ū)
# x_prev = deepcopy(x̄)
x_prev = deepcopy(xref) # TODO deepcopy(x̄)

# ## problem
prob = IterativeLQR.problem_data(model, obj, cons)

for ρ in [1e-3]#, 3e-4]#, 1e-4]#, 3e-5]#, 1e-5, 3e-6, 1e-6]
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
	    max_iter=25,
		max_al_iter=3,
	    ρ_init=1.0,
	    ρ_scale=10.0)

	# Update initialization trajectory
	x_prev, u_prev = IterativeLQR.get_trajectory(prob)
	visualize(env, x_prev)
end

scatter(abs.(goal(prob.m_data.x[T], zeros(0), zeros(0))))
# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)
@show norm(vcat([ctrl_lmt(prob.m_data.x[t], prob.m_data.u[t], zeros(0)) for t=1:T-1]...), Inf)

# jldsave(joinpath(@__DIR__, "rexhopper_traj.jld2"), x_sol=x_sol, u_sol=u_sol)

plot([x[5] for x in x_sol])


# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(env.vis, cam_pos=[0,-3,2], zoom=3)
open(env.vis)

set_floor!(env.vis, x=6.0, y=6.0, z=0.01, alt=0.0, color=RGBA(0.5,0.5,0.5,1.0))
set_camera!(env.vis, cam_pos=[5,-6,6], zoom=3)
set_light!(env.vis)

convert_frames_to_video_and_gif("atlas_tryna_jump")

render_static(env.vis)
# open(joinpath(@__DIR__, "atlas_6_steps.html"), "w") do file
#     write(file, static_html(env.vis))
# end
