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
global const ROTATE = Ref{Float64}(0.0)
loss(:box2d, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[loss(:box2d, pairs0, data0 + [i;zeros(19)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[loss(:box2d, pairs0, data0 + [0;i;zeros(18)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 4 contact locations [y,z] -> 9 params in total
################################################################################
using Optim

function termination_callback(opt; value_tol=1e-4)
	return opt.value < value_tol
end

solver = LBFGS(;m=100,
        alphaguess=Optim.LineSearches.InitialStatic(),
        linesearch=Optim.LineSearches.HagerZhang(),
        P = 1e2*I(9),
		)

lower = [0.00, +0.05, +0.05, +0.05, -1.00, -1.00, +0.05, -1.00, -1.00]
upper = [0.80, +1.00, +1.00, +1.00, -0.05, -0.05, +1.00, -0.05, -0.05]

function d2data(d)
	cf = d[1]
	data = [cf; 0.05; 0; +d[2]; +d[3];
			cf; 0.05; 0; +d[4]; +d[5];
			cf; 0.05; 0; +d[6]; +d[7];
			cf; 0.05; 0; +d[8]; +d[9];
			]
	return data
end

∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(9))

function fg!(F,G,d)
	# do common computations here
	l, ∇ = loss(:box2d, pairs0, d2data(d), n_sample=50)
	if G != nothing
		G .= ∇d2data' * ∇
	end
	if F != nothing
		value = l
	    return value
	end
end
d0 = [0.40, +0.50, +0.50, +0.50, -0.50, -0.50, +0.50, -0.50, -0.50]
ROTATE[] = 0.0
optimize(Optim.only_fg!(fg!), lower, upper, d0, Fminbox(solver),
	Optim.Options(
		callback = termination_callback,
		show_trace = true),
	; inplace = false,
	)


# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a single side length.
################################################################################
using Optim

solver = LBFGS(;m=100,
        alphaguess=Optim.LineSearches.InitialStatic(),
        linesearch=Optim.LineSearches.HagerZhang(),
        P = 1e2*I(2),
		)

lower = [0.00, 0.05]
upper = [0.80, 1.00]

function d2data(d)
	cf, side = d
	data = [cf; 0.05; 0; +side; +side;
			cf; 0.05; 0; +side; -side;
			cf; 0.05; 0; -side; +side;
			cf; 0.05; 0; -side; -side;
			]
	return data
end

∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

function fg!(F,G,d)
	# do common computations here
	l, ∇ = loss(:box2d, pairs0, d2data(d), n_sample=50)
	if G != nothing
		G .= ∇d2data' * ∇
	end
	if F != nothing
		value = l
	    return value
	end
end

d0 = [0.40, 0.50]
ROTATE[] = 0.0
optimize(Optim.only_fg!(fg!), lower, upper, d0, Fminbox(solver),
	Optim.Options(
		callback = termination_callback,
		show_trace = true),
	; inplace = false,
	)



# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

################################################################################
# Visualization
################################################################################
