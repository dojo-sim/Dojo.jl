# Load packages
using Dojo
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include( "../utils.jl")

mech = getmechanism(:box2d, Δt=0.05, g=-9.81, cf=0.2, radius = 0.05, side = 0.50);
initialize!(mech, :box2d, x=[-1,1.], v=[2,1.], θ=0.1, ω=2.)
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact = true)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0.2], [1,0.4]],
				   :vlims => [-2ones(2), ones(2)],
				   :θlims => [-π, π],
				   :ωlims => [-10, 10])
mech_kwargs = Dict(:cf => 0.1, :radius => 0.05, :side => 0.5)
generate_dataset(:box2d, H=0.75, N=15,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs)


################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:box2d; N = 15, mech_kwargs...)

data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
clean_loss(:box2d, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[clean_loss(:box2d, pairs0, data0 + [i;zeros(27)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:box2d, pairs0, data0 + [0;i;zeros(26)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 4 contact locations [y,z] -> 9 params in total
################################################################################


include("../quasi_newton.jl")
function d2data(d)
	cf = d[1]
	data = [cf; 0;0; 0.05; 0; +d[2]; +d[3];
			cf; 0;0; 0.05; 0; +d[4]; +d[5];
			cf; 0;0; 0.05; 0; +d[6]; +d[7];
			cf; 0;0; 0.05; 0; +d[8]; +d[9];
			]
	return data
end

∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(9))

d0 = [0.40, +0.50, +0.50, +0.50, -0.50, -0.50, +0.50, -0.50, -0.50]
lower = [0.00, +0.05, +0.05, +0.05, -1.00, -1.00, +0.05, -1.00, -1.00]
upper = [0.80, +1.00, +1.00, +1.00, -0.05, -0.05, +1.00, -0.05, -0.05]

function f0(d; rot=0)
	return clean_loss(:box2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
# dsol = bfgs_solve(f0, g0, H0, d0, iter=200, gtol=1e-8, ftol=3e-5,
dsol = quasi_newton_solve(f0, fgH0, d0, iter=200, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-6)

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
	cf, side = d
	data = [cf; 0;0; 0.05; 0; +side; +side;
			cf; 0;0; 0.05; 0; +side; -side;
			cf; 0;0; 0.05; 0; -side; +side;
			cf; 0;0; 0.05; 0; -side; -side;
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 0.75]
lower = [0.00, 0.05]
upper = [0.80, 1.00]

function f0(d; rot=0)
	return clean_loss(:box2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box2d, pairs0, d2data(d), n_sample=15, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
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
