using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
Pkg.instantiate()

# ## Setup
using Dojo
using Plots
using Random
using MeshCat
using DojoEnvironments

# ## Open visualizer
vis = Visualizer()
render(vis)

# ## Include methods
methods_dir = joinpath("../../methods")
include(joinpath(methods_dir, "filename.jl"))
include(joinpath(methods_dir, "initial_state.jl"))
include(joinpath(methods_dir, "quasi_newton.jl"))
include(joinpath(methods_dir, "dataset.jl"))
include(joinpath(methods_dir, "loss.jl"))


# ## ---------------------------------------------------------------------------
# parameters
# ## ---------------------------------------------------------------------------
S = 1
# we multiplied all the distance by length_scaling for better data scaling
# we must rescale the gravity term accordingly since gravity contains length (g == kg.m.s^-2)
length_scaling = 20.0
timestep = 1/148 * S
gravity = -9.81 * scaling
friction_coefficient = 0.16
radius = 0.00
side = 2.00
mode = :box
opts_step = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
opts_grad = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
model = :block
N = 100

mech_kwargs = Dict(
	:timestep => timestep,
	:gravity => gravity,
	:friction_coefficient => friction_coefficient,
	:radius => radius,
	:side => side,
	:mode => mode)

# ## ---------------------------------------------------------------------------
# simulation example
# ## ---------------------------------------------------------------------------
mech = get_mechanism(model; mech_kwargs...);
initialize!(mech, model,
	position=[0.0, -1.0, 1.0],
	velocity=[0.0, 2.0, 1.0],
	angular_velocity=[2.0, 5.0, 10.0])
storage = simulate!(mech, 5.0, record=true,
    opts=opts_step)
visualize(mech, storage, vis=vis, show_contact = true)

# ## ---------------------------------------------------------------------------
# Generate & Save Dataset
# ## ---------------------------------------------------------------------------
generate_hardware_dataset(model,
	N=N,
	opts=opts_step,
	mech_kwargs=mech_kwargs,
	vis=vis,
	sleep_ratio=0.01,
	S=S,
	)

# ## ---------------------------------------------------------------------------
# Load Dataset
# ## ---------------------------------------------------------------------------
params0, trajs0 = open_dataset(model; experiment_type="hardware", N=N, mech_kwargs...)
data0 = params0[:data]
nd = sum(data_dim.(mech.contacts))
data_contacts0 = data0[end-nd+1:end]


# ## ---------------------------------------------------------------------------
# rollout
# ## ---------------------------------------------------------------------------
data_storage = trajs0[1]
z0 = get_maximal_state(data_storage, 1)

mech_kwargs = Dict(
	:timestep => timestep,
	:gravity => gravity,
	:friction_coefficient => friction_coefficient,
	:radius => radius,
	:side => side,
	:mode => mode)


mech = get_mechanism(model; mech_kwargs...)
set_maximal_state!(mech, z0)
storage = simulate!(mech, length(data_storage) * timestep, opts=opts_step)
vis, anim = visualize(mech, storage, vis=vis, name=:robot, color=RGBA(0.3, 0.3, 0.3, 1.0))
vis, anim = visualize(mech, data_storage, vis=vis, animation=anim, name=:data)


# ## ---------------------------------------------------------------------------
# Optimization Objective: Evaluation & Gradient
# ## ---------------------------------------------------------------------------
indices0 = 70:70
function f0(d; rot=0, n_sample=0, trajs=trajs0, N=N, indices=indices0)
	f = 0.0
	mechanism = get_mechanism(model; mech_kwargs...)
	for i = 1:N
		fi, Z = loss(mechanism, d_to_data_contacts(d), trajs[i], indices,
			opts_step=opts_step, opts_grad=opts_grad, derivatives=false)
		f += fi
		# visualize(mechanism, generate_storage(mechanism, Z), vis=vis, animation=anim)
		# sleep(1.0)
	end
	return f
end

function fgH0(d; rot=0, n_sample=0, trajs=trajs0, N=N, indices=indices0)
	mechanism = get_mechanism(model; mech_kwargs...)
	f = 0.0
	nd = sum(data_dim.(mechanism.contacts))
	g = zeros(nd)
	H = zeros(nd,nd)
	for i = 1:N
		fi, gi, Hi = loss(mech, d_to_data_contacts(d), trajs[i], indices,
			opts_step=opts_step, opts_grad=opts_grad, derivatives=true)
		f += fi
		g += gi
		H += Hi
	end
	return f, data_mask' * g, data_mask' * H * data_mask
end


# ## ---------------------------------------------------------------------------
# We learn a single coefficient of friction and a single side length -> 2 params in total
# ## ---------------------------------------------------------------------------
function d_to_data_contacts(d)
	friction_coefficient = d[1]
	half_side = d[2]
	data_contacts = [
			friction_coefficient; 0.0; +half_side; +half_side; -half_side;
			friction_coefficient; 0.0; +half_side; -half_side; -half_side;
			friction_coefficient; 0.0; -half_side; +half_side; -half_side;
			friction_coefficient; 0.0; -half_side; -half_side; -half_side;
			friction_coefficient; 0.0; +half_side; +half_side; +half_side;
			friction_coefficient; 0.0; +half_side; -half_side; +half_side;
			friction_coefficient; 0.0; -half_side; +half_side; +half_side;
			friction_coefficient; 0.0; -half_side; -half_side; +half_side;
			]
	return data_contacts
end

data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))

# ## ---------------------------------------------------------------------------
# cost Landscape
# ## ---------------------------------------------------------------------------
d1 = [0.17, 0.9735]
d_to_data_contacts(d1) - data_contacts0
f0(d1)
f1, g1, H1 = fgH0(d1)
- H1 \ g1

X = 0.0:0.04:0.8
F1 = [f0([x; d1[2:end]]) for x in X]
plot(X, F1)

X = 1.6:0.02:2.00
F2 = [f0([d1[1]; x; d1[3:end]]) for x in X]
plot(X, F2)


# ## ---------------------------------------------------------------------------
# solve
# ## ---------------------------------------------------------------------------
d0 = [0.40, 2.0]
lower = [0.00, 0.05]
upper = [0.80, 2.00]

indices0 = 70:70
dsol = quasi_newton_solve(f0, fgH0, d0, iter=20, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end



# ## ---------------------------------------------------------------------------
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
# ## ---------------------------------------------------------------------------
function d_to_data_contacts(d)
	friction_coefficient = d[1]
	data_contacts = [
			friction_coefficient; 0; +d[2:4];
			friction_coefficient; 0; +d[5:7];
			friction_coefficient; 0; +d[8:10];
			friction_coefficient; 0; +d[11:13];
			friction_coefficient; 0; +d[14:16];
			friction_coefficient; 0; +d[17:19];
			friction_coefficient; 0; +d[20:22];
			friction_coefficient; 0; +d[23:25];
			]
	return data_contacts
end

data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(25))

# ## ---------------------------------------------------------------------------
# cost Landscape
# ## ---------------------------------------------------------------------------
d1 = [0.17,
	+0.95, +0.95, -0.95,
	+0.95, -0.95, -0.95,
	-0.95, +0.95, -0.95,
	-0.95, -0.95, -0.95,
	+0.95, +0.95, +0.95,
	+0.95, -0.95, +0.95,
	-0.95, +0.95, +0.95,
	-0.95, -0.95, +0.95,
	]
f0(d1)
f1, g1, H1 = fgH0(d1)
- H1 \ g1

X = 0.0:0.05:0.8
F1 = [f0([x; d1[2:end]]) for x in X]
plot(X, F1)

X = 0.7:0.05:1.1
F2 = [f0([d1[1]; x; d1[3:end]]) for x in X]
plot(X, F2)


# ## ---------------------------------------------------------------------------
# solve
# ## ---------------------------------------------------------------------------
d0 = [0.40,
	+2.00, +2.00, -2.00,
	+2.00, -2.00, -2.00,
	-2.00, +2.00, -2.00,
	-2.00, -2.00, -2.00,
	+2.00, +2.00, +2.00,
	+2.00, -2.00, +2.00,
	-2.00, +2.00, +2.00,
	-2.00, -2.00, +2.00,
	]
lower = [0.00,
	+0.05, +0.05, -2.00,
	+0.05, -2.00, -2.00,
	-2.00, +0.05, -2.00,
	-2.00, -2.00, -2.00,
	+0.05, +0.05, +0.05,
	+0.05, -2.00, +0.05,
	-2.00, +0.05, +0.05,
	-2.00, -2.00, +0.05,
	]
upper = [0.80,
	+2.00, +2.00, -0.05,
	+2.00, -0.05, -0.05,
	-0.05, +2.00, -0.05,
	-0.05, -0.05, -0.05,
	+2.00, +2.00, +2.00,
	+2.00, -0.05, +2.00,
	-0.05, +2.00, +2.00,
	-0.05, -0.05, +2.00,
	]

indices0 = 50:52
dsol = quasi_newton_solve(f0, fgH0, d0, iter=20, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)


losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end

# ## ---------------------------------------------------------------------------
# save solution for visualization purposes
# ## ---------------------------------------------------------------------------
jldsave(joinpath(@__DIR__, "../..", "data", "hardware", "visualization",
	"hardware_block_solution.jld2"), dsol=dsol)
