# using Pkg
# Pkg.activate(joinpath(@__DIR__, "../../.."))
# Pkg.instantiate()

# ## Setup
# Toss data from https://github.com/DAIRLab/contact-nets, distances scaled by factor 20
using Dojo
using DojoEnvironments
using Plots
using JLD2
using ForwardDiff
using LinearAlgebra

# ## Include methods
include("utils.jl")


# ## Parameters
# we multiplied all the distance by distance_scaling for better data scaling
# we must rescale the gravity term accordingly since gravity contains length (g == kg.m.s^-2)
distance_scaling = 20.0
opts = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
model = :block
N = 100

mech_kwargs = Dict(
	:timestep => 1/148, :gravity => -9.81 * distance_scaling,
	:friction_coefficient => 0.16, :edge_length => 0.1 * distance_scaling)

# ## Load dataset
dataset = jldopen(joinpath(@__DIR__, "data", "datasets", "real_block.jld2"))
storages = dataset["storages"]
JLD2.close(dataset)

function parameter_stack(θ)
	# [friction_coefficient; contact_radius; contact_origin]
	return [
		θ[1]; 0; +θ[2:4];
		θ[1]; 0; +θ[5:7];
		θ[1]; 0; +θ[8:10];
		θ[1]; 0; +θ[11:13];
		θ[1]; 0; +θ[14:16];
		θ[1]; 0; +θ[17:19];
		θ[1]; 0; +θ[20:22];
		θ[1]; 0; +θ[23:25];
	]
end

data_mask = ForwardDiff.jacobian(x -> parameter_stack(x), zeros(25))

# ## Optimization Objective: Evaluation & Gradient
timesteps = 50:52
function f0(d; storages=storages, timesteps=timesteps)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	for i = 1:100
		fi = loss(mechanism, parameter_stack(d), storages[i], timesteps;
			opts, derivatives=false)
		f += fi
	end
	return f
end

function fgH0(d; storages=storages, timesteps=timesteps)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	nd = sum(Dojo.data_dim.(mechanism.contacts))
	g = zeros(nd)
	H = zeros(nd,nd)
	for i = 1:100
		fi, gi, Hi = loss(mechanism, parameter_stack(d), storages[i], timesteps;
			opts, derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, data_mask' * g, data_mask' * H * data_mask
end

# ## Initial guess
guess = [0.40,
	+2.00, +2.00, -2.00,
	+2.00, -2.00, -2.00,
	-2.00, +2.00, -2.00,
	-2.00, -2.00, -2.00,
	+2.00, +2.00, +2.00,
	+2.00, -2.00, +2.00,
	-2.00, +2.00, +2.00,
	-2.00, -2.00, +2.00,
]
guess_error = f0(guess)

# ## Solve
lower_bound = [0.00,
	+0.05, +0.05, -2.00,
	+0.05, -2.00, -2.00,
	-2.00, +0.05, -2.00,
	-2.00, -2.00, -2.00,
	+0.05, +0.05, +0.05,
	+0.05, -2.00, +0.05,
	-2.00, +0.05, +0.05,
	-2.00, -2.00, +0.05,
]
upper_bound = [0.80,
	+2.00, +2.00, -0.05,
	+2.00, -0.05, -0.05,
	-0.05, +2.00, -0.05,
	-0.05, -0.05, -0.05,
	+2.00, +2.00, +2.00,
	+2.00, -0.05, +2.00,
	-0.05, +2.00, +2.00,
	-0.05, -0.05, +2.00,
]

solution = quasi_newton_solve(f0, fgH0, guess; iter=20, gtol=1e-8, ftol=1e-6,
	lower_bound, upper_bound, reg=1e-9)

# ## Result
solution_parameters = solution[1]
solution_error = f0(solution_parameters)

# ## Visualize
vis = Visualizer()
storage = storages[1]

mech_kwargs = Dict(
	:timestep => 1/148, :gravity => -9.81 * distance_scaling, :color => RGBA(1,0,0,0.5),
	:friction_coefficient => 0.16, :edge_length => 0.1 * distance_scaling)
mech = get_mechanism(model; mech_kwargs...)
vis, animation = visualize(mech, storage; vis, return_animation=true, name=:original)


position = storage.x[1][1]-[0;0;solution_parameters[2]]
orientation = storage.q[1][1]
velocity = storage.v[1][1]
angular_velocity = storage.ω[1][1]
edge_length = 2*sum(abs.(solution_parameters[2:25]))/24 # mean edge length for visualization

mech_kwargs = Dict(
	:timestep => 1/148, :gravity => -9.81 * distance_scaling, :color => RGBA(0,1,0,0.5),
	:friction_coefficient => solution_parameters[1], :edge_length => edge_length)
mech = get_mechanism(model; mech_kwargs...)
set_data!(mech.contacts, parameter_stack(solution_parameters)) # set actual learned contact points

initialize!(mech, model; position, velocity, orientation, angular_velocity)
storage = simulate!(mech, length(storage)/148; record=true, opts)
visualize(mech, storage; vis, animation, name=:learned)
