# Load packages
using Dojo
using Plots
using Random
using MeshCat
using OSFLoader

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include( "../methods/utils.jl")
include( "../methods/quasi_newton.jl")
include( "bunny_utils.jl")

aprilcube_options = ColliderOptions(
	impact_damper=3e5,
	impact_spring=3e4,
	sliding_drag=0.00,
	sliding_friction=0.1,
	rolling_drag=0.05,
	rolling_friction=0.01,
	coulomb_smoothing=3e1,
	coulomb_regularizer=1e-3)
mech = get_mechanism(:aprilcube, timestep=1/148, gravity=-9.81*10, side=1.0, options=aprilcube_options);
get_data(mech.contacts)
initialize!(mech, :aprilcube, position=[0,-1,1.], velocity=3.0*[0,2,1.], angular_velocity=4.0*[2,2,2.])
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=3e-4, rtol=3e-4, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)


################################################################################
# Generate & Save Dataset
################################################################################
generate_hardware_dataset(N=10, sleep_ratio=0.1, scale=0.5)


################################################################################
# Load Dataset
################################################################################
params0, trajs0 = open_dataset(:aprilcube; N=10)
data0 = params0[:data]
data_contacts0 = data0[end-4:end]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
timestep = params0[:timestep]
gravity = params0[:g]
model = :aprilcube
indices0 = 30:80

function f0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity, options=aprilcube_options)
	for i = 1:N
		f += loss(mechanism, d_to_data_contacts(d), trajs[i], indices,
			opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=5, indices=indices0)
	mechanism = get_mechanism(model, timestep=timestep, gravity=gravity, options=aprilcube_options)
	f = 0.0
	g = zeros(5)
	H = zeros(5,5)
	for i = 1:N
		fi, gi, Hi = loss(mech, d_to_data_contacts(d), trajs[i], indices,
			opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, data_mask' * g, data_mask' * H * data_mask
end

mechanism = get_mechanism(model, timestep=timestep, gravity=gravity, options=deepcopy(aprilcube_options))
f = loss(mechanism, d_to_data_contacts([0.1]), storage, indices0,
	opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)



################################################################################
# Optimization Algorithm: Quasi Newton:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################
# F = [f0([x]) for x in 0:0.04:0.5]
plot(0:0.04:0.5, F)
# G = [fgH0([x])[2][1] for x in 0:0.04:0.5]
plot(0:0.04:0.5, G)
# f0([0.0])
# fgH0([0.0])

z0 = get_maximal_state(storage, 50)
get_contact_gradients!(mech, z0, [2.30, 0.1, 0,0,0])
_, _, J1 = get_contact_gradients!(mech, z0, [2.30, 0.1, 0,0,0])
J1


function d_to_data_contacts(d)
	bounciness = data_contacts0[1]
	friction_coefficient = d[1]
	# bounciness = d[1]
	# friction_coefficient = d[2]
	data_contacts = [bounciness; friction_coefficient; data_contacts0[3:5]]
	return data_contacts
end
data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(1))
# data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))

# d0 = [2.0, 0.0]
# lower = [-3.0, 0.0]
# upper = [+3.0, 1.0]

d0 = [0.4]
lower = [0.0]
upper = [1.0]

# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=10, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9, α0=1.0)


losses = f0.(dsol[2])
plt = plot()
plot!(plt, losses, ylims=[0,10], xlabel="iterations", ylabel="prediction loss", legend=false)
scatter!(plt, losses, ylims=[0,10], xlabel="iterations", ylabel="prediction loss", legend=false)

for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end

################################################################################
# Visualization
################################################################################
id = 2
position_init = trajs0[id].x[1][1]
velocity_init = trajs0[id].v[1][1]
orientation_init = trajs0[id].q[1][1]
angular_velocity_init = trajs0[id].ω[1][1]
z_init = [position_init; velocity_init; vector(orientation_init); angular_velocity_init]

mech = get_mechanism(:aprilcube, timestep=timestep, gravity=gravity, options=aprilcube_options);
# set_data!(mech.contacts, [dsol[2][1]; data_contacts0[3:5]])
set_data!(mech.contacts, [data_contacts0[1]; dsol[2][1]; data_contacts0[3:5]])
set_maximal_state!(mech, z_init)
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=3e-4, rtol=3e-4, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, color=RGBA(1,1,1,1.), name=:initial)

aprilcube_options = ColliderOptions(
	impact_damper=3e5,
	impact_spring=3e4,
	sliding_drag=0.00,
	sliding_friction=0.1,
	rolling_drag=0.05,
	rolling_friction=0.01,
	coulomb_smoothing=3e1,
	coulomb_regularizer=1e-3)
mech = get_mechanism(:aprilcube, timestep=timestep, gravity=gravity, options=aprilcube_options);
# set_data!(mech.contacts, [dsol[1]; data_contacts0[3:5]])
set_data!(mech.contacts, [data_contacts0[1]; dsol[1]; data_contacts0[3:5]])

set_maximal_state!(mech, z_init)
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=3e-4, rtol=3e-4, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.7,0.7,0.7,1.), name=:learned)

mech = get_mechanism(:aprilcube, timestep=timestep, gravity=gravity, friction_coefficient=0.4);
vis, anim = visualize(mech, trajs0[id], vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,1.), name=:robot)

storage_init = generate_storage(mech, [z_init])
vis, anim = visualize(mech, storage_init, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,0.3), name=:start)


# convert_frames_to_video_and_gif("aprilcube_learning_friction")

################################################################################
# DEBUG
################################################################################
