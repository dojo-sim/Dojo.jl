# Load packages
using Dojo
using Plots
using Random
using MeshCat
using DojoEnvironments

# Open visualizer
vis = Visualizer()
render(vis)

# Include new files
methods_dir = joinpath("../../methods")
include(joinpath(methods_dir, "filename.jl"))
include(joinpath(methods_dir, "initial_state.jl"))
include(joinpath(methods_dir, "data.jl"))
include(joinpath(methods_dir, "data_jacobian.jl"))
include(joinpath(methods_dir, "quasi_newton.jl"))
include(joinpath(methods_dir, "dataset.jl"))
include(joinpath(methods_dir, "loss.jl"))

# parameters
timestep = 0.05
gravity = -9.81
radius = 0.05
side = 0.50
friction_coefficient = 0.1
opts_step = SolverOptions(btol=3e-4, rtol=3e-4)
opts_grad = SolverOptions(btol=3e-4, rtol=3e-4)
model = :block2d
N = 15

mech = get_mechanism(model,
	timestep=timestep,
	gravity=gravity,
	friction_coefficient=friction_coefficient,
	radius=radius,
	side=side);
initialize!(mech, model,
	position=[-1,1.], linear_velocity=[2,1.],
	orientation=0.1,
	angular_velocity=2.)
storage = simulate!(mech, 5.0, record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact = true)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0.2], [1,0.4]],
				   :vlims => [-2ones(2), ones(2)],
				   :θlims => [-π, π],
				   :ωlims => [-10, 10])
mech_kwargs = Dict(
	:friction_coefficient => friction_coefficient,
	:radius => radius,
	:side => side)
generate_dataset(model,
	H=0.75,
	N=N,
	opts=opts_step,
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	vis=vis,
	sleep_ratio=0.01,
	)


################################################################################
# Load Dataset
################################################################################
params0, trajs0 = open_dataset(model; N=N, mech_kwargs...)
data0 = params0[:data]
data_contacts0 = data0[end-4*5+1:end]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
indices0 = 1:N-1
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=10, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model; mech_kwargs...)
	for i = 1:N
		fi, Z = loss(mechanism, d_to_data_contacts(d), trajs[i], indices,
			step_opts=step_opts, opts_grad=opts_grad, derivatives=false)
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
			step_opts=step_opts, opts_grad=opts_grad, derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, data_mask' * g, data_mask' * H * data_mask
end

################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 4 contact locations [y,z] -> 9 params in total
################################################################################
function d_to_data_contacts(d)
	friction_coefficient = d[1]
	data_contacts = [friction_coefficient; 0.05; 0; +d[2]; +d[3];
			friction_coefficient; 0.05; 0; +d[4]; +d[5];
			friction_coefficient; 0.05; 0; +d[6]; +d[7];
			friction_coefficient; 0.05; 0; +d[8]; +d[9];
			]
	return data_contacts
end


data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(9))

X = 0.0:0.02:0.5
F = [f0([x; 0.25; 0.25; 0.25; -0.25; -0.25; 0.25; -0.25; -0.25]) for x in X]
plot(X, F)

d0 = [0.40, +0.50, +0.50, +0.50, -0.50, -0.50, +0.50, -0.50, -0.50]
lower = [0.00, +0.05, +0.05, +0.05, -1.00, -1.00, +0.05, -1.00, -1.00]
upper = [0.80, +1.00, +1.00, +1.00, -0.05, -0.05, +1.00, -0.05, -0.05]

# Main.@profiler
dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-6)

losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end


d0
f0(d0)
f0(dsol)

# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a single side length.
################################################################################
include("../quasi_newton.jl")
function d2data(d)
	friction_coefficient, side = d
	data = [friction_coefficient; 0;0; 0.05; 0; +side; +side;
			friction_coefficient; 0;0; 0.05; 0; +side; -side;
			friction_coefficient; 0;0; 0.05; 0; -side; +side;
			friction_coefficient; 0;0; 0.05; 0; -side; -side;
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 0.75]
lower = [0.00, 0.05]
upper = [0.80, 1.00]

function f0(d; rot=0)
	return clean_loss(:block2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:block2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

d0
f0(d0)
f0(dsol)
f0([0.1, 0.5])

# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

################################################################################
# Visualization
################################################################################
