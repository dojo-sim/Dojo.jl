# Load packages
using Dojo
using Plots
using Random
using MeshCat
using OSFLoader
# using ForwardDiff

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
methods_dir = joinpath("../../system_identification/methods")
include(joinpath(methods_dir, "filename.jl"))
include(joinpath(methods_dir, "initial_state.jl"))
include(joinpath(methods_dir, "data.jl"))
include(joinpath(methods_dir, "data_jacobian.jl"))
include(joinpath(methods_dir, "quasi_newton.jl"))
include("methods/dataset.jl")
include("methods/loss.jl")


mech = get_mechanism(:nerf, nerf=:halfsoap, timestep=0.01, gravity=-9.81, friction_coefficient=0.05);
mech.contacts[1].model.collision.collider.options =
	ColliderOptions(
	impact_damper=1e5,
	impact_spring=3e4,
	sliding_drag=0.1,
	sliding_friction=0.05,
	rolling_drag=0.0,
	rolling_friction=0.01,
	coulomb_smoothing=3e1,
	coulomb_regularizer=1e-3,)

initialize!(mech, :nerf,
	position=[0,-1,0.5],
	velocity=[0,5,1.],
	orientation=Quaternion(-0.2, 0.7, 0.6, -0.2, false),
	angular_velocity=[0,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:nerf => :halfsoap,
				   :xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [[-3,-3,-0.5], [3,3,0.]],
				   :ωlims => [-1ones(3), 1ones(3)])
mech_kwargs = Dict(:nerf => :halfsoap,
				   :friction_coefficient => friction_coefficient)

generate_dataset(:nerf,
	N=10,
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
params0, trajs0 = open_dataset(:nerf; N=50, mech_kwargs...)
data0 = params0[:data]
data_contacts0 = data0[end-4:end]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
timestep = 0.01
gravity = -9.81
model = :nerf
nerf = :halfsoap
indices0 = 80:90
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model, nerf=nerf, timestep=timestep, gravity=gravity)
	for i = 1:N
		f += loss(mechanism, d_to_data_contacts(d), trajs[i], indices, opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	mechanism = get_mechanism(model, nerf=nerf, timestep=timestep, gravity=gravity)
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
function d_to_data_contacts(d)
	bounciness = data_contacts0[1]
	friction_coefficient = d[1]
	data_contacts = [bounciness; friction_coefficient; data_contacts0[3:5]]
	# bounciness = d[1]
	# friction_coefficient = d[2]
	# data_contacts = [bounciness; friction_coefficient; data_contacts0[3:5]]
	return data_contacts
end
# data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))
data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(1))

F = [f0([x]) for x in 0:0.02:1]
plot(0:0.02:1, F)


# d0 = [-2.30, 0.10]
# lower = [-3.0, 0.0]
# upper = [+3.0, 1.0]

d0 = [0.8]
lower = [0.0]
upper = [1.0]

# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=1000, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end
################################################################################
# Visualization
################################################################################
mech = get_mechanism(:nerf, nerf=nerf, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, [data_contacts0[1]; dsol[2][1]; data_contacts0[3:5]])
initialize!(mech, :nerf,
	position=[0,-1,0.5],
	velocity=[0,5,1.],
	orientation=Quaternion(-0.2, 0.7, 0.6, -0.2, false),
	angular_velocity=[0,2,2.])
storage = simulate!(mech, 3.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, color=RGBA(1,1,1,1.), name=:initial)


mech = get_mechanism(:nerf, nerf=nerf, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, [data_contacts0[1]; dsol[1]; data_contacts0[3:5]])
# mech.bodies[1].mass /= 10
# mech.bodies[1].inertia /= 10
initialize!(mech, :nerf,
	position=[0,-1,0.5],
	velocity=[0,5,1.],
	orientation=Quaternion(-0.2, 0.7, 0.6, -0.2, false),
	angular_velocity=[0,2,2.])
storage = simulate!(mech, 3.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.7,0.7,0.7,1.), name=:learned)

mech = get_mechanism(:nerf, nerf=nerf, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, data_contacts0)
initialize!(mech, :nerf,
	position=[0,-1,0.5],
	velocity=[0,5,1.],
	orientation=Quaternion(-0.2, 0.7, 0.6, -0.2, false),
	angular_velocity=[0,2,2.])
storage = simulate!(mech, 3.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,1.), name=:robot)

z_init = get_maximal_state(storage, 1)
storage_init = generate_storage(mech, [z_init])
vis, anim = visualize(mech, storage_init, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,0.3), name=:start)


convert_frames_to_video_and_gif("halfsoap_learning_friction")
