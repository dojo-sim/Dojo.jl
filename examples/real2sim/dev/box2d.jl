# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))
include(joinpath(module_dir(), "examples", "real2sim", "utils.jl"))

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
generate_dataset(:box2d, H=1.0, N=15,
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
loss(mech, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[loss(mech, pairs0, [0.1, 0.5+i, 0,0,0], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[loss(mech, pairs0, [0.1+i, 0.5, 0,0,0], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: L-BFGS
################################################################################

# Solution for bad cost landscape
	# use L-BFGS
	# use longer horizons to compute the loss (currently the horizon is 1 step)
		# maybe we just need to sum the gradients along the horizon
		# maybe we need to chain them together using th chain rule (not sure how stable this is)

using Optim
solver = LBFGS(; m = 100,
        alphaguess = Optim.LineSearches.InitialStatic(),
        linesearch = Optim.LineSearches.HagerZhang(),
        P = 1e1*I(2),
        precondprep = (P, x) -> nothing,
        manifold = Optim.Flat(),
		)

lower = [0.0, 0.0]
upper = [0.80, 2.0]
f(d) = loss(mech, pairs0, [d;zeros(3)])[1]
g(d) = loss(mech, pairs0, [d;zeros(3)])[2][1:2]
d0 = [0.40, 1.0]
optimize(f, g, lower, upper, d0, Fminbox(solver),
	Optim.Options(x_tol = 1e-4,
		f_tol = 1e-7,
		g_tol = 1e-7,
		show_trace = true),
	; inplace = false,
	)

# We can learn the coefficient of friction and the radius form 15*0.75 seconds of
# recording. We use the simulator to evaluate the loss and its gradients by differentiating
# through the simulator. With gradient information we can use L-BFGS

################################################################################
# Visualization
################################################################################
