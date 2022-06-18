using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
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
include(joinpath(module_dir(), "DojoEnvironments", "src", "sphere", "deps", "texture.jl"))

################################################################################
# parameters
################################################################################
timestep = 0.02
gravity = -9.81
friction_coefficient = 0.1
radius = 0.50
opts_step = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
opts_grad = SolverOptions(btol=3e-4, rtol=3e-4, undercut=3.0)
model = :sphere
N = 1

mech_kwargs = Dict(
	:timestep => timestep,
	:gravity => gravity,
	:friction_coefficient => friction_coefficient,
	:radius => radius)


################################################################################
# simulation example
################################################################################
mech = get_mechanism(model; mech_kwargs...);
initialize!(mech, model,
	position=[0,-1,2.],
	velocity=[0,2,1.],
	angular_velocity=[2,5,10.],
	)
storage = simulate!(mech, 5.0, record=true,
    opts=opts_step)
visualize(mech, storage, vis=vis)
sphere_texture!(vis, mech)



################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0], [1,1,0.2]],
				   :vlims => [-1ones(3), 1ones(3)],
				   :ωlims => [-5ones(3), 5ones(3)])

generate_dataset(model,
	H=0.70,
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
nd = sum(data_dim.(mech.contacts))
data_contacts0 = data0[end-nd+1:end]


################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
indices0 = 10:12
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


################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 4 contact locations [y,z] -> 9 params in total
################################################################################
function d_to_data_contacts(d)
	friction_coefficient = d[1]
	contact_radius = d[2]
	data_contacts = [friction_coefficient; contact_radius; 0;0;0]
	return data_contacts
end

data_mask = FiniteDiff.finite_difference_jacobian(d -> d_to_data_contacts(d), zeros(2))

d1 = [0.10, 0.5]
f0(d1)
f1, g1, H1 = fgH0(d1)
- H1 \ g1

X = 0.0:0.05:0.8
F1 = [f0([x; d1[2:end]]) for x in X]
plot(X, F1)

X = 0.0:0.05:0.55
F2 = [f0([d1[1]; x; d1[3:end]]) for x in X]
plot(X, F2)



d0 = [0.00, +1.00]
lower = [0.00, +0.05]
upper = [0.80, +1.00]

dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)



losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end


# We can learn the coefficient of friction and the side dimenson of the cube
# form 100*0.70 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS
