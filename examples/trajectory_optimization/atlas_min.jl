# PREAMBLE

# PKG_SETUP

# ## setup
using Dojo
using IterativeLQR
using LinearAlgebra 

# ## visualizer
vis = Visualizer() 
open(vis)

# ## system
include(joinpath(@__DIR__, "../../environments/atlas/methods/template.jl"))

gravity = -9.81
timestep = 0.05
friction_coefficient = 0.8
damper = 50.0
spring = 0.0
model_type = :armless
env = get_environment(:atlas,
    representation=:minimal,
    timestep=timestep,
    gravity=gravity,
    friction_coefficient=friction_coefficient,
    damper=damper,
    spring=spring,
	model_type=model_type)

# ## dimensions
n = env.num_states
m = env.num_inputs
d = 0

# ## reference trajectory
N = 1
initialize!(env.mechanism, :atlas)
xref = atlas_trajectory(env.mechanism; 
	timestep=timestep, 
	r=0.0, 
	x=0.03, 
	z=0.85, 
	N=12, 
	Ncycles=N)
zref = [minimal_to_maximal(env.mechanism, x) for x in xref]
visualize(env, xref)

## gravity compensation TODO: solve optimization problem instead
mech = get_mechanism(:atlas, 
	timestep=timestep, 
	gravity=gravity, 
	friction_coefficient=friction_coefficient, 
	damper=1000.0,
	spring=spring, 
	model_type=model_type)

initialize!(mech, :atlas)
storage = simulate!(mech, 0.10, 
	record=true, 
	verbose=false)

visualize(mech, storage, 
	vis=env.vis)

F_damper = get_damper_impulses(env.mechanism)
u_damper = F_damper * env.mechanism.timestep
u_control = u_damper[6 .+ (1:15)]

mech = get_mechanism(:atlas, 
	timestep=timestep, 
	gravity=gravity, 
	friction_coefficient=friction_coefficient, 
	damper=0.0,
	spring=spring, 
	model_type=model_type)

function controller!(mechanism, k)
    set_input!(mechanism, u_damper)
    return
end

initialize!(mech, :atlas, 
	body_position=[0,0,0.0])
storage = simulate!(mech, 0.50, controller!, 
	record=true, 
	verbose=false)

visualize(mech, storage, 
	vis=env.vis)

function get_damper_impulses(mechanism::Mechanism{T}) where T
	joints = mechanism.joints
	# set the controls in the equality constraints
	off = 0
	nu = input_dimension(mechanism)
	u = zeros(nu)
	for joint in joints
		pbody = get_body(mechanism, joint.parent_id)
		if typeof(pbody) <: Body
			F = damper_impulses(mechanism, joint, pbody)
			oF = 0
			for joint in [joint.translational, joint.rotational]
				nf, nF = size(nullspace_mask(joint))
				u[off .+ (1:nf)] .= nullspace_mask(joint) * F[oF .+ (1:nF)]
				off += nf
				oF += nF
			end
		else
			for joint in [joint.translational, joint.rotational]
				nf, nF = size(nullspace_mask(joint))
				off += nf
			end
		end
	end
	return u
end

# ## horizon
T = N * (25 - 1) + 1

# ## model
dyn = IterativeLQR.Dynamics(
    (y, x, u, w) -> dynamics(y, env, x, u, w),
    (dx, x, u, w) -> dynamics_jacobian_state(dx, env, x, u, w),
    (du, x, u, w) -> dynamics_jacobian_input(du, env, x, u, w),
    n, n, m, d)

model = [dyn for t = 1:T-1]

# ## rollout
x1 = xref[1]
ū = [u_control for t = 1:T-1]
w = [zeros(d) for t = 1:T-1]
x̄ = IterativeLQR.rollout(model, x1, ū, w)
visualize(env, x̄)

# ## objective
qt = [0.3; 0.05; 0.05; 0.01 * ones(3); 0.01 * ones(3); 0.01 * ones(3); fill(0.002, 30)...]
ots = [(x, u, w) -> transpose(x - xref[t]) * Diagonal(timestep * qt) * (x - xref[t]) + transpose(u) * Diagonal(timestep * 0.002 * ones(m)) * u for t = 1:T-1]
oT = (x, u, w) -> transpose(x - xref[end]) * Diagonal(timestep * qt) * (x - xref[end])

cts = IterativeLQR.Cost.(ots, n, m, d)
cT = IterativeLQR.Cost(oT, n, 0, 0)
obj = [cts..., cT]

# ## constraints
function goal(x, u, w)
    Δ = x - xref[end]
	Δ[3] -= 0.40
    return Δ[collect(1:3)]
end

cont = IterativeLQR.Constraint()
conT = IterativeLQR.Constraint(goal, n, 0)
cons = [[cont for t = 1:T-1]..., conT]

# ## problem
prob = IterativeLQR.problem_data(model, obj, cons)
IterativeLQR.initialize_controls!(prob, ū)
IterativeLQR.initialize_states!(prob, x̄)

# ## solve
IterativeLQR.solve!(prob,
    verbose=true,
	linesearch=:armijo,
    α_min=1.0e-5,
    obj_tol=1.0e-3,
    grad_tol=1.0e-3,
    max_iter=100,
    max_al_iter=5,
    ρ_init=1.0,
    ρ_scale=10.0)

vis= Visualizer()
open(env.vis)

# ## solution
x_sol, u_sol = IterativeLQR.get_trajectory(prob)
@show IterativeLQR.eval_obj(prob.m_data.obj.costs, prob.m_data.x, prob.m_data.u, prob.m_data.w)
@show prob.s_data.iter[1]
@show norm(goal(prob.m_data.x[T], zeros(0), zeros(0)), Inf)

# ## visualize
x_view = [[x_sol[1] for t = 1:15]..., x_sol..., [x_sol[end] for t = 1:15]...]
visualize(env, x_view)

set_camera!(env.vis, 
	zoom=5.0, 
	cam_pos=[0.0, -5.0, 0.0])

z = [minimal_to_maximal(env.mechanism, x) for x in x_sol]
t = 1 #10, 20, 30, 41
set_robot(env.vis, env.mechanism, z[t])
