# Load packages
using Dojo
using Plots
using Random
using MeshCat
using ForwardDiff

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include( "../methods/utils.jl")
include( "../methods/quasi_newton.jl")
include( "bunny_utils.jl")

mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.3);

initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)
θ0 = get_data(mech.contacts)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [[-4,-4,-0.5], [4,4,0.]],
				   :ωlims => [-2ones(3), 2ones(3)])
mech_kwargs = Dict(:friction_coefficient => 0.3)

generate_dataset(:bunny,
	N = 5,
	opts=SolverOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	show_contact=false,
	sleep_ratio=1.0,
	vis=vis,
	)

################################################################################
# Load Dataset
################################################################################
params0, trajs0 = open_dataset(:bunny; N=5, mech_kwargs...)
data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
timestep = 0.01
gravity = -9.81
model = :bunny
indices0 = 20:199
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity)
	for i = 1:N
		f += loss(mechanism, d2data(d), trajs[i], indices, opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity)
	f = 0.0
	g = zeros(4)
	H = zeros(4,4)
	for i = 1:N
		fi, gi, Hi = loss(mech, d2data(d), trajs[i], indices, opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end


################################################################################
# Optimization Algorithm: Quasi Newton:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################


function d2data(d)
	friction_coefficient = d[1]
	collider_origin = d[2:4]
	data = [friction_coefficient; collider_origin]
	return data
end
∇d2data = FiniteDiff.finite_difference_jacobian(d -> d2data(d), zeros(4))

d0 = [0.05,0,0,0.0]
lower = [0.0, -0.2, -0.2, -0.2]
upper = [1.0, +0.2, +0.2, +0.2]


# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=10, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

dsol
θ0


################################################################################
# Visualization
################################################################################
mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, dsol[2][1])
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, color=RGBA(1,1,1,1.), name=:initial)

mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, dsol[1])
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.7,0.7,0.7,1.), name=:learned)

mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, θ0)
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,1.), name=:bunny)


# convert_frames_to_video_and_gif("bunny_learning_friction_com")

################################################################################
# DEBUG
################################################################################

xa = sones(3)
qa = Quaternion(1,0,0,0.0)
xb = sones(3)
qb = Quaternion(1,0,0,0.0)
η = szeros(1)

joint = mech.joints[1].translational
relative = :parent
impulse_map(relative, joint, xa, qa, xb, qb, η)
@benchmark $impulse_map($relative, $joint, $xa, $qa, $xb, $qb, $η)

joint = mech.joints[1].rotational
relative = :parent
impulse_map(relative, joint, xa, qa, xb, qb, η)
@benchmark $impulse_map($relative, $joint, $xa, $qa, $xb, $qb, $η)

joint = mech.joints[1].translational
relative = :child
impulse_map(relative, joint, xa, qa, xb, qb, η)
@benchmark $impulse_map($relative, $joint, $xa, $qa, $xb, $qb, $η)

joint = mech.joints[1].rotational
relative = :child
impulse_map(relative, joint, xa, qa, xb, qb, η)
@benchmark $impulse_map($relative, $joint, $xa, $qa, $xb, $qb, $η)



function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    relative = (body.id == contact.parent_id ? :parent : :child)
    pbody = get_body(mechanism, contact.parent_id)
    cbody = get_body(mechanism, contact.child_id)
    return impulse_map(relative, contact.model, pbody, cbody, mechanism.timestep)
end


contact = mech.contacts[1]
body = mech.bodies[1]
impulse_map(mech, contact, body)
using BenchmarkTools
@benchmark $impulse_map($mech, $contact, $body)

relative = (body.id == contact.parent_id ? :parent : :child)
pbody = get_body(mech, contact.parent_id)
cbody = get_body(mech, contact.child_id)


timestep = 0.01
model = contact.model
relative=:parent
impulse_map(relative, model, pbody, cbody, timestep)
@benchmark $impulse_map($relative, $model, $pbody, $cbody, $timestep)

relative=:child
impulse_map(relative, model, pbody, cbody, timestep)
@benchmark $impulse_map($relative, $model, $pbody, $cbody, $timestep)
