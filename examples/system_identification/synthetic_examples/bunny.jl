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
mech.contacts[1].model.collision.collider.options =
	ColliderOptions(
	impact_damper=3e5,
	impact_spring=3e4,
	sliding_drag=0.1,
	sliding_friction=0.3,
	rolling_drag=0.0,
	rolling_friction=0.01,
	coulomb_smoothing=3e1,
	coulomb_regularizer=1e-3,)

initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [[-4,-4,-0.5], [4,4,0.]],
				   :ωlims => [-2ones(3), 2ones(3)])
mech_kwargs = Dict(:friction_coefficient => 0.3)

generate_dataset(:bunny,
	N = 50,
	opts=SolverOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	show_contact=false,
	sleep_ratio=0.0,
	vis=vis,
	)

################################################################################
# Load Dataset
################################################################################
params0, trajs0 = open_dataset(:bunny; N=50, mech_kwargs...)
data0 = params0[:data]
data_contacts0 = data0[end-4:end]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
timestep = 0.01
gravity = -9.81
model = :bunny
indices0 = 80:90
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity)
	for i = 1:N
		f += loss(mechanism, d_to_data_contacts(d), trajs[i], indices, opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity)
	f = 0.0
	g = zeros(5)
	H = zeros(5,5)
	for i = 1:N
		fi, gi, Hi = loss(mech, d_to_data_contacts(d), trajs[i], indices, opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, data_mask' * g, data_mask' * H * data_mask
end


################################################################################
# Optimization Algorithm: Quasi Newton:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################
F = [f0([2.30, x]) for x in 0:0.02:1]
plot(0:0.02:1, F)

function d_to_data_contacts(d)
	bounciness = d[1]
	friction_coefficient = d[2]
	data_contacts = [bounciness; friction_coefficient; data_contacts0[3:5]]
	return data_contacts
end
data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))

d0 = [-2.30, 0.10]
lower = [-3.0, 0.0]
upper = [+3.0, 1.0]

# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=1000, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

################################################################################
# Visualization
################################################################################
mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, [dsol[2][1]; data_contacts0[3:5]])
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, color=RGBA(1,1,1,1.), name=:initial)

mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, [dsol[1]; data_contacts0[3:5]])
mech.bodies[1].mass /= 10
mech.bodies[1].inertia /= 10
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.7,0.7,0.7,1.), name=:learned)

mech = get_mechanism(:bunny, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, θ0)
initialize!(mech, :bunny, position=[0,-1,1.], velocity=[0,2,1.], angular_velocity=[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,1.), name=:robot)


# convert_frames_to_video_and_gif("bunny_learning_friction_com")

################################################################################
# DEBUG
################################################################################
