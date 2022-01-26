# Load packages
using Dojo
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "env", "sphere", "deps", "texture.jl"))
include( "../utils.jl")

mech = getmechanism(:sphere, timestep=0.05, g=-9.81, radius=0.5, cf=0.1);
initialize!(mech, :sphere, x=[0,0,0.3], v=[0,0.5,0.], ω=[10,0,0.])
storage = simulate!(mech, 0.5, record=true, verbose=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6))
visualize(mech, storage, vis=vis)
sphere_texture!(vis, mech)

################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0], [1,1,0.2]],
				   :vlims => [-1ones(3), 1ones(3)],
				   :ωlims => [-5ones(3), 5ones(3)])
mech_kwargs = Dict(:cf => 0.1, :radius => 0.5)
generate_dataset(:sphere, H=0.75, N=15,
	opts=SolverOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs)


################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:sphere; N=15, mech_kwargs...)

data0 = params0[:data]

################################################################################
# CLEAN Optimization Objective: Evaluation & Gradient
################################################################################

clean_loss(:sphere, pairs0, data0, opts=SolverOptions(btol=3e-4, rtol=3e-4))

[clean_loss(:sphere, pairs0, [0.1, 0,0, 0.5+i, 0,0,0], opts=SolverOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:sphere, pairs0, [0.1+i, 0,0, 0.5, 0,0,0], opts=SolverOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# CLEAN Optimization Algorithm: L-BFGS
################################################################################

# Solution for bad cost landscape
	# use L-BFGS
	# use longer horizons to compute the loss (currently the horizon is 1 step)
		# maybe we just need to sum the gradients along the horizon
		# maybe we need to chain them together using th chain rule (not sure how stable this is)

include("../quasi_newton.jl")
function d2data(d)
	data = [d[1]; 0; 0; d[2]; 0; 0; 0]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 1.0]
lower = [0.0, 0.0]
upper = [0.80, 2.0]

function f0(d; rot=0)
	return clean_loss(:sphere, pairs0, d2data(d), n_sample=15, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:sphere, pairs0, d2data(d), n_sample=15, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end

dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=3e-5,
	H0=1e1*Matrix(Diagonal(ones(2))), lower=lower, upper=upper)

d0
f0(d0)
f0(dsol)
f0([0.1, 0.50])

using Plots; pyplot()
x=range(0,stop=0.80,length=10)
y=range(0,stop=0.75,length=10)
f(x,y) = log(10, clean_loss(:sphere, pairs0, d2data([x,y]),
	opts=SolverOptions(btol=3e-4, rtol=3e-4))[1])
plot(x,y,f,st=:surface,camera=(20,40))

# We can learn the coefficient of friction and the radius form 15*0.75 seconds of
# recording. We use the simulator to evaluate the loss and its gradients by differentiating
# through the simulator. With gradient information we can use L-BFGS

################################################################################
# Visualization
################################################################################
