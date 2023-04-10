# ### Setup
# PKG_SETUP
using Dojo
using Plots
using DojoEnvironments
using JLD2
using LinearAlgebra

# ### Include methods
include("utilities.jl")

# ### Parameters
model = :sphere
N = 10
timesteps = 10:12 # at which trajectory will be used for learning

mech_kwargs = Dict(
	:timestep => 0.02, :gravity => -9.81,
	:friction_coefficient => 0.2, :radius => 0.50)

# ### Generate dataset
mech = get_mechanism(model; mech_kwargs...)

xlims = [[0,0,0], [1,1,0.2]]
vlims = [-1*ones(3), 1*ones(3)]
ωlims = [-5*ones(3), 5*ones(3)]

storages = []

for i = 1:N
	position = xlims[1] + rand(3) .* (xlims[2] - xlims[1])
	velocity = vlims[1] + rand(3) .* (vlims[2] - vlims[1])
	angular_velocity = ωlims[1] + rand(3) .* (ωlims[2] - ωlims[1])
	initialize!(mech, model; position, velocity, angular_velocity)
	storage = simulate!(mech, 2; record=true)
	push!(storages, storage)
	## visualize(mech, storage)
end


# ### Save dataset
jldsave(joinpath(@__DIR__, "data", "datasets", "synthetic_sphere.jld2"); storages)

# ### Save and load dataset
dataset = jldopen(joinpath(@__DIR__, "data", "datasets", "synthetic_sphere.jld2"))
storages = dataset["storages"]
JLD2.close(dataset)

# ### Optimization Objective: Evaluation & Gradient
function parameter_stack(θ)
	## [friction_coefficient; contact_radius; contact_origin]
	return [θ; zeros(3)]
end

function f0(θ; storages=storages, timesteps=timesteps)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	for i = 1:length(storages)
		fi = loss(mechanism, parameter_stack(θ), storages[i], timesteps; derivatives=false)
		f += fi
	end
	return f
end

function fgH0(θ; storages=storages, timesteps=timesteps)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	nd = sum(Dojo.data_dim.(mechanism.contacts))
	g = zeros(nd)
	H = zeros(nd,nd)

	for i = 1:length(storages)
		fi, gi, Hi = loss(mech, parameter_stack(θ), storages[i], timesteps; derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, g[1:2], H[1:2,1:2] # We only care about the first two parameters
end

# ### Cost Landscape
X = 0:0.02:0.4
Y = 0.47:0.005:0.53
surface(X, Y, (X, Y) -> f0([X, Y]; storages))


# ### Initial guess
guess = [0.0, 1.0]
guess_error = f0(guess)

# ### Solve
lower_bound = [0, 0.05]
upper_bound = [0.8, 1]
solution = quasi_newton_solve(f0, fgH0, guess; iter=20, gtol=1e-8, ftol=1e-6,
	lower_bound, upper_bound, reg=1e-9)

# ### Result
solution_parameters = solution[1]
solution_error = f0(solution_parameters)

# ### Visualize
vis = Visualizer()
storage = storages[1]

mech_kwargs = Dict(
	:timestep => 0.02, :gravity => -9.81, :color => RGBA(1,0,0,0.5),
	:friction_coefficient => 0.2, :radius => 0.50)
mech = get_mechanism(model; mech_kwargs...)
vis, animation = visualize(mech, storage; vis, return_animation=true, name=:original)

position = storage.x[1][1]-[0;0;solution_parameters[2]]
orientation = storage.q[1][1]
velocity = storage.v[1][1]
angular_velocity = storage.ω[1][1]

mech_kwargs = Dict(
	:timestep => 0.02, :gravity => -9.81, :color => RGBA(0,1,0,0.5),
	:friction_coefficient => solution_parameters[1], :radius => solution_parameters[2])
mech = get_mechanism(model; mech_kwargs...)

initialize!(mech, model; position, velocity, angular_velocity)
storage = simulate!(mech, 2; record=true)
visualize(mech, storage; vis, animation, name=:learned)
render(vis)