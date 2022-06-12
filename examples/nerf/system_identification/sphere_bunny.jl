# using Pkg
# Pkg.develop(path=joinpath(@__DIR__, "../../../OSFLoader"))

# Load packages
using Dojo
using Plots
using Random
using MeshCat
# using ForwardDiff

# Open visualizer
vis = Visualizer()
render(vis)

# Include new files
methods_dir = joinpath("../../system_identification/methods")
include(joinpath(methods_dir, "filename.jl"))
include(joinpath(methods_dir, "initial_state.jl"))
include(joinpath(methods_dir, "data.jl"))
include(joinpath(methods_dir, "data_jacobian.jl"))
include(joinpath(methods_dir, "quasi_newton.jl"))
# include(joinpath(methods_dir, "dataset.jl"))
# include(joinpath(methods_dir, "loss.jl"))
include("methods/dataset.jl")
include("methods/loss_body_contact.jl")

################################################################################
# Parameters
################################################################################
timestep = 0.01
gravity = -9.81
torsional_friction = 0.01
sliding_friction = 0.20
sphere_mass = 10.0
sphere_radius = 0.25
r_sphere = 1.0
v_sphere = 4.0
nerf_position = [0,0,0.35]
collider_options = ColliderOptions(
	impact_damper=3e5,
	impact_spring=3e4,
	sliding_drag=0.01,
	sliding_friction=sliding_friction,
	torsional_drag=0.00,
	torsional_friction=torsional_friction,
	rolling_drag=0.00,
	rolling_friction=0.01,
	coulomb_smoothing=3e1,
	coulomb_regularizer=1e-3,)

################################################################################
# Simulation
################################################################################
mech = get_mechanism(:nerf_sphere, nerf=:bunny,
	mass = sphere_mass,
	radius = sphere_radius,
	timestep=timestep,
	gravity=gravity,
	friction_coefficient=sliding_friction,
	collider_options=collider_options);
mech.contacts[1]
q_nerf = Quaternion(normalize(rand(4))...)
α_sphere = 2π * rand()

initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [cos(α_sphere), sin(α_sphere),0.],
	sphere_velocity=v_sphere * [-cos(α_sphere), -sin(α_sphere),0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:v_sphere => v_sphere,
				   :r_sphere => r_sphere,
				   :nerf_position => nerf_position,
				   )
mech_kwargs = Dict(:nerf => :bunny,
				   :mass => sphere_mass,
	               :radius => sphere_radius,
				   :timestep => timestep,
				   :gravity => gravity,
				   :friction_coefficient => sliding_friction,
				   :collider_options => collider_options,
				   )

generate_dataset(:nerf_sphere,
	N=50,
	opts=SolverOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	show_contact=false,
	sleep_ratio=0.8,
	vis=vis,
	)

################################################################################
# Load Dataset
################################################################################
params0, trajs0 = open_dataset(:nerf_sphere; N=50, mech_kwargs...)
data0 = params0[:data]

data_body_contact0 = [
	data0[17:23]; # body 1
	data0[37:43]; # body 2
	data0[end-14:end]] # contact 1, 2, 3

get_data(mech.bodies[1])
get_data(mech.bodies[2])

get_data(mech.contacts[1])
get_data(mech.contacts[2])
get_data(mech.contacts[3])

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
model = :nerf_sphere
indices0 = 70:100
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=15, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model; mech_kwargs...)
	for i = 1:N
		fi, _ = loss(mechanism, d_to_data(d), trajs[i], indices,
			opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
		f += fi
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=15, indices=indices0)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	g = zeros(29)
	H = zeros(29,29)
	for i = 1:N
		fi, gi, Hi = loss(mech, d_to_data(d), trajs[i], indices,
			opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=true)
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
function d_to_data_body_contact(d)
	data_body_contact = deepcopy(data_body_contact0)

	data_body_contact[1] = d[1] # mass
	sliding_friction = d[2]
	data_body_contact[14 .+ [2,12]] .= sliding_friction

	# data_body_contact[1] = d[1] # mass

	return data_body_contact
end

function d_to_data(d)
	data = deepcopy(data0)
	data_body_contact = d_to_data_body_contact(d)

	data[16 .+ (1:7)] .= data_body_contact[1:7]
	data[16 + 20 .+ (1:7)] .= data_body_contact[7 .+ (1:7)]
	data[end-14:end] .= data_body_contact[end-14:end]
	return data
end

d_to_data_body_contact(rand(2))

data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_body_contact(d), zeros(2))

F = [f0([0.25, x]) for x in 0:0.02:0.5]
plot(0:0.02:0.5, F)
F = [f0([x, 0.00]) for x in 1:0.2:5]
plot(1:0.2:5, F)

# F = [f0([x]) for x in 0.25:0.15:10]
# plot(0.25:0.15:10, F)

d0 = [0.25, 0.0]
lower = [0.25, 0.0]
upper = [10.0, 0.5]

# d0 = [0.25]
# lower = [0.25]
# upper = [10.0]

# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=20, gtol=1e-8, ftol=1e-6,
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
initialize!(mech, :nerf, position=[0,-1,1.], velocity=[0,3,3.], angular_velocity=[1,2,1.])
storage = simulate!(mech, 6.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, color=RGBA(1,1,1,1.), name=:initial)

mech = get_mechanism(:nerf, nerf=nerf, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, [data_contacts0[1]; dsol[1]; data_contacts0[3:5]])
# mech.bodies[1].mass /= 10
# mech.bodies[1].inertia /= 10
initialize!(mech, :nerf, position=[0,-1,1.], velocity=[0,3,3.], angular_velocity=[1,2,1.])
storage = simulate!(mech, 6.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.7,0.7,0.7,1.), name=:learned)

mech = get_mechanism(:nerf, nerf=nerf, timestep=0.01, gravity=-9.81, friction_coefficient=0.4);
set_data!(mech.contacts, data_contacts0)
initialize!(mech, :nerf, position=[0,-1,1.], velocity=[0,3,3.], angular_velocity=[1,2,1.])
storage = simulate!(mech, 6.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,1.), name=:robot)

z_init = get_maximal_state(storage, 1)
storage_init = generate_storage(mech, [z_init])
vis, anim = visualize(mech, storage_init, vis=vis, animation=anim, color=RGBA(0.2,0.2,0.2,0.3), name=:start)

# render_static(vis)
# open("/home/simon/bunny_system_identification.html", "w") do file
#     write(file, static_html(vis))
# end



# convert_frames_to_video_and_gif("bunny_learning_friction")
