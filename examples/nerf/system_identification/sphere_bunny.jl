# using Pkg
# Pkg.develop(path=joinpath(@__DIR__, "../../../OSFLoader"))

# Load packages
using Dojo
using Plots
using Random
using MeshCat
using DojoEnvironments
# using ForwardDiff

# Open visualizer
vis = Visualizer()
open(vis)
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
sliding_friction = 0.10
sphere_mass = 10.0
sphere_radius = 0.25
r_sphere = 1.0
v_sphere = 2.0
nerf_position = [0,0.20,0.35]
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

function ctrl!(mechanism::Mechanism{T}, k::Int; kp=1*2e-0, kd=1*3e-1,
		goal_position=[0,0,0.], goal_orientation=[0,0,0.]) where T

	timestep = mechanism.timestep

	sphere = get_body(mechanism, :sphere)
    x_sphere = current_position(sphere.state)
	axis, angle = axis_angle(current_orientation(sphere.state))
	θ_sphere = angle .* axis
    v_sphere, ω_sphere = current_velocity(sphere.state)

	u_position = (goal_position - x_sphere) * kp - kd * v_sphere
    u_orientation = (goal_orientation - θ_sphere) * kp - kd * ω_sphere
	input = [szeros(6); u_position; u_orientation] / timestep
	set_input!(mechanism, input)
	return input
end

q_nerf = Quaternion(normalize(rand(4))...)
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=[0,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)

################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:r_sphere => r_sphere,
				   :v_sphere => v_sphere,
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

# generate_dataset(:nerf_sphere,
# 	N=50,
# 	ctrl! = ctrl!,
# 	opts=SolverOptions(btol=3e-4, rtol=3e-4),
# 	init_kwargs=init_kwargs,
# 	mech_kwargs=mech_kwargs,
# 	show_contact=false,
# 	sleep_ratio=0.006,
# 	vis=vis,
# 	)

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
indices0 = 50:60
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=10, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model; mech_kwargs...)
	for i = 1:N
		fi, Z = loss(mechanism, d_to_data(d), trajs[i], indices,
			opts=SolverOptions(btol=3e-4, rtol=3e-4), derivatives=false)
		f += fi
		# visualize(mechanism, generate_storage(mechanism, Z), vis=vis, animation=anim)
		# sleep(1.0)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=10, indices=indices0)
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

	mass0 = data_body_contact0[1]
	inertia0 = data_body_contact0[2:7]

	mass = d[1]
	data_body_contact[1] = mass # mass
	data_body_contact[2:7] .= mass / mass0 * inertia0 # inertia
	sliding_friction = d[2]
	data_body_contact[14 .+ [2,12]] .= sliding_friction

	# mass = d[1]
	# data_body_contact[1] = mass # mass
	# sliding_friction = d[2]
	# data_body_contact[14 .+ [2,12]] .= sliding_friction

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

vis, anim = visualize(mech, trajs0[7], vis=vis, name=:ref)
f0([data0[17], 0.10])
fgH0([data0[17], 0.10])

F1 = [f0([3.22, x]) for x in 0:0.05:0.4]
plot(0:0.05:0.4, F1)
F1 = [f0([3.00, x]) for x in 0:0.05:0.4]
plot(0:0.05:0.4, F1)
F1 = [f0([0.25, x]) for x in 0:0.05:0.4]
plot(0:0.05:0.4, F1)
F2 = [f0([x, 0.10]) for x in 1.2:0.4:5]
plot(1.2:0.4:5, F2)
F2 = [f0([x, 0.00]) for x in 1.2:0.4:5]
plot(1.2:0.4:5, F2)

# F = [f0([x]) for x in 0.25:0.15:10]
# plot(0.25:0.15:10, F)

d0 = [0.25, 0.01]
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
render(vis)
scale_vis = 0.25
q_nerf = Quaternion(normalize([1,0.5,0,0])...)
offset = [-0.60,0,0]

# initial
mech = get_mechanism(:nerf_sphere, nerf=:bunny; mech_kwargs...)
for body in mech.bodies
	body.shape.scale *= scale_vis
end
set_data!(mech, d_to_data(dsol[2][1]))
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=v_sphere * [-1,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
initial_storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
for x in initial_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end
vis, anim = visualize(mech, initial_storage, vis=vis,
	color=RGBA(1,1,1,1.), name=:initial)

# learned
mech = get_mechanism(:nerf_sphere, nerf=:bunny; mech_kwargs...)
for body in mech.bodies
	body.shape.scale *= scale_vis
end
set_data!(mech, d_to_data(dsol[1]))
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=v_sphere * [-1,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
learned_storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
for x in learned_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end
vis, anim = visualize(mech, learned_storage, vis=vis, animation=anim,
	color=RGBA(0.7,0.7,0.7,1.), name=:learned)

# ground_truth
mech = get_mechanism(:nerf_sphere, nerf=:bunny; mech_kwargs...)
for body in mech.bodies
	body.shape.scale *= scale_vis
end
set_data!(mech, data0)
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=v_sphere * [-1,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
ground_truth_storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
for x in ground_truth_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end
vis, anim = visualize(mech, ground_truth_storage, vis=vis, animation=anim,
	color=RGBA(0.2,0.2,0.2,1.), name=:ground_truth)

################################################################################
# Export
################################################################################
# vis = Visualizer()
render(vis)
panda_mech = get_mechanism(:panda, damper=100.0, model_type=:end_effector, contact=true,
	joint_limits=[[-10.0, -1.7628, -2.8973, -0.0698, -2.8973, -3.7525, -2.8973, -0.00],
				  [ 10.0,  1.7628,  2.8973,  3.0718,  2.8973,  0.0175,  2.8973,  0.04]],)

initialize!(panda_mech, :panda, joint_angles=[π, -0.8, 0.0, 1.6, 0.0, -3.2, 0.0, 0.0, 0.0])
storage = simulate!(panda_mech, 0.01)
visualize(panda_mech, storage, vis=vis)
q_end_effector = current_orientation(get_body(panda_mech, 12).state)

# visual_storage = deepcopy(initial_storage)
# visual_storage = deepcopy(learned_storage)
visual_storage = deepcopy(ground_truth_storage)

z_panda = panda_inverse_kinematics_trajectory(panda_mech,
	visual_storage.x[2][1:200],
	[q_end_effector for i=1:200])

panda_storage = generate_storage(panda_mech, z_panda)
vis, anim = visualize(panda_mech, panda_storage, vis=vis, animation=anim, name=:panda,
	show_contact=true)


α = 1.0
colors = [RGBA(0.2,0.2,0.2,α), RGBA(0.6,0.6,0.6,α), RGBA(1.0,1.0,1.0,α)];
index = [1,100,200]
i = 3

panda_storage_i = generate_storage(panda_mech, [get_maximal_state(panda_storage, index[i]) for t=1:200])
vis, anim = visualize(panda_mech, panda_storage_i, vis=vis, animation=anim, name=:panda_i,
	show_contact=true, color=colors[i])

visual_storage_i = generate_storage(mech, [get_maximal_state(visual_storage, index[i]) for t=1:200])
vis, anim = visualize(mech, visual_storage_i, vis=vis, animation=anim, name=:sphere_bunny_i,
	color=colors[i])

# set_camera!(vis, zoom=5)


# jldsave(joinpath(@__DIR__, "sphere_bunny.jld2"), x=storage.x[2], q=[vector(q) for q in storage.q[2]])
# jldopen(joinpath(@__DIR__, "sphere_bunny.jld2"))["x"]
# jldopen(joinpath(@__DIR__, "sphere_bunny.jld2"))["q"]


################################################################################
# Dataset Visualization
################################################################################
vis = Visualizer()
open(vis)

params0, trajs0 = open_dataset(:nerf_sphere; N=50, mech_kwargs...)

scale_vis = 0.25
offset = [-0.60,0,0]

data_storage = trajs0[21]
for x in data_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end

z_panda = panda_inverse_kinematics_trajectory(panda_mech,
	data_storage.x[2][1:200],
	[q_end_effector for i=1:200])

panda_storage = generate_storage(panda_mech, z_panda)
vis, anim = visualize(panda_mech, panda_storage, vis=vis, name=:panda,
	show_contact=true)

vis, anim = visualize(mech, data_storage, vis=vis, animation=anim, name=:sphere_bunny,
	color=RGBA(0.9, 0.9, 0.9, 1))

for i = 1:16
	convert_frames_to_video_and_gif("sphere_bunny_dataset_$i")
end


################################################################################
# Learning Visualization
################################################################################
vis = Visualizer()
render(vis)
open(vis)

scale_vis = 0.25
offset = [-0.60,0,0]

mech = get_mechanism(:nerf_sphere, nerf=:bunny; mech_kwargs...)
for body in mech.bodies
	body.shape.scale *= scale_vis
end

# learned
set_data!(mech, d_to_data(dsol[2][end]))
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=v_sphere * [-1,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
learned_storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
for x in learned_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end
vis, anim = visualize(mech, learned_storage, vis=vis,
	color=RGBA(0.9,0.9,0.9,1.), name=:learned)

# ground_truth
set_data!(mech, data0)
initialize!(mech, :nerf_sphere,
	nerf_position=nerf_position,
	sphere_position=r_sphere * [1,0,0.],
	sphere_velocity=v_sphere * [-1,0,0.],
	nerf_orientation=q_nerf,
	sphere_orientation=one(Quaternion),
	)
ground_truth_storage = simulate!(mech, 2.0, ctrl!, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
for x in ground_truth_storage.x
	x .*= scale_vis
	x .+= fill(offset, 200)
end
vis, anim = visualize(mech, ground_truth_storage, vis=vis, animation=anim,
	color=RGBA(0.2,0.2,0.2,1.), name=:ground_truth)


z_learned_panda = panda_inverse_kinematics_trajectory(panda_mech,
	learned_storage.x[2][1:200],
	[q_end_effector for i=1:200])
z_ground_truth_panda = panda_inverse_kinematics_trajectory(panda_mech,
	ground_truth_storage.x[2][1:200],
	[q_end_effector for i=1:200])

panda_learned_storage = generate_storage(panda_mech, z_learned_panda)
panda_ground_truth_storage = generate_storage(panda_mech, z_ground_truth_panda)

vis, anim = visualize(panda_mech, panda_learned_storage, vis=vis, animation=anim, name=:learned,
	show_contact=true, color=RGBA(0.9,0.9,0.9,1.))
vis, anim = visualize(panda_mech, panda_ground_truth_storage, vis=vis, animation=anim, name=:ground_truth,
	show_contact=true, color=RGBA(0.2,0.2,0.2,1.))

settransform!(vis[:learned], Translation(0,0.3,0.0))
set_floor!(vis, x=3.5, y=3.5)

convert_frames_to_video_and_gif("sphere_bunny_interactive_perception_20")
