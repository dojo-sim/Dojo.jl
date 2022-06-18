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


################################################################################
# parameters
################################################################################
timestep = 0.02
gravity = -9.81
friction_coefficient = 0.1
radius = 0.00
side = 0.50
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

################################################################################
# simulation example
################################################################################
mech = get_mechanism(model; mech_kwargs...);
initialize!(mech, model,
	position=[0,-1,1.],
	velocity=[0,2,1.],
	angular_velocity=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=opts_step)
visualize(mech, storage, vis=vis, show_contact = true)

################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [-2ones(3), [2,2,-1.]],
				   :ωlims => [-6ones(3), 6ones(3)])

generate_dataset(model,
	H=0.40,
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
indices0 = 10:10
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

d1 = [0.11,
	+0.25, +0.25, -0.25,
	+0.25, -0.25, -0.25,
	-0.25, +0.25, -0.25,
	-0.25, -0.25, -0.25,
	+0.25, +0.25, +0.25,
	+0.25, -0.25, +0.25,
	-0.25, +0.25, +0.25,
	-0.25, -0.25, +0.25,
	]
f0(d1)
f1, g1, H1 = fgH0(d1)
- H1 \ g1

X = 0.0:0.05:0.8
F1 = [f0([x; d1[2:end]]) for x in X]
plot(X, F1)

X = 0.0:0.05:0.30
F2 = [f0([d1[1]; x; d1[3:end]]) for x in X]
plot(X, F2)



d0 = [0.40,
	+0.5, +0.5, -0.5,
	+0.5, -0.5, -0.5,
	-0.5, +0.5, -0.5,
	-0.5, -0.5, -0.5,
	+0.5, +0.5, +0.5,
	+0.5, -0.5, +0.5,
	-0.5, +0.5, +0.5,
	-0.5, -0.5, +0.5,
	]
lower = [0.00,
	+0.05, +0.05, -1.00,
	+0.05, -1.00, -1.00,
	-1.00, +0.05, -1.00,
	-1.00, -1.00, -1.00,
	+0.05, +0.05, +0.05,
	+0.05, -1.00, +0.05,
	-1.00, +0.05, +0.05,
	-1.00, -1.00, +0.05,
	]
upper = [0.80,
	+1.00, +1.00, -0.05,
	+1.00, -0.05, -0.05,
	-0.05, +1.00, -0.05,
	-0.05, -0.05, -0.05,
	+1.00, +1.00, +1.00,
	+1.00, -0.05, +1.00,
	-0.05, +1.00, +1.00,
	-0.05, -0.05, +1.00,
	]


dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)


losses = f0.(dsol[2])
for (i,l) in enumerate(losses)
	println("($(i-1),$(l/losses[1]))")
end


# We can learn the coefficient of friction and the side dimenson of the cube
# form 100*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator.

################################################################################
# Visualization
################################################################################

# Visualize Dataset with solution
datasol = deepcopy(data0)
datasol[end-nd+1:end] = d_to_data_contacts(dsol[1])
set_data!(mech, datasol)
for traj in trajs0[1:11]
	visualize(mech, traj, vis=vis, show_contact = true)
	sleep(0.5)
end

# Visualize Progress made during 'learning'
for (i,d) in enumerate(dsol[2][1:8])
	datasol = deepcopy(data0)
	datasol[end-nd+1:end] = d_to_data_contacts(d)
	set_data!(mech, datasol)
	visualize(mech, trajs0[i], vis=vis, show_contact = true)
	sleep(1.5)
end
