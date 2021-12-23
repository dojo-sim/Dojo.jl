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

mech = getmechanism(:box, Δt=0.05, g=-9.81, cf=0.2, radius=0.00, side=0.50, mode=:box);
initialize!(mech, :box, x=[0,-1,1.], v=[0,2,1.], ω=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact = false)


#
# ineqc1 = mech.ineqconstraints.values[1]
# ∂g∂simdata(mech, ineqc1)
#
# simmat, simmatFD = simdata_jacobian(mech)
# plot(Gray.(abs.(simmat)))
# plot(Gray.(abs.(simmatFD)))
# plot(Gray.(abs.(1e110 .* simmatFD)))
#
# norm((simmat + simmatFD)[4:7, 1:7], Inf)
# norm((simmat + simmatFD)[7:14, 1:7], Inf)
#
# length(mech.ineqconstraints.values[1])
# length(mech.bodies.values[1])





################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [-2ones(3), [2,2,-1.]],
				   :ωlims => [-6ones(3), 6ones(3)])
mech_kwargs = Dict(:cf => 0.1, :radius => 0.0, :side => 0.5)
generate_dataset(:box, H=0.40, N=25,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	sleep_ratio=0.01,
	show_contact=false)


################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:box; N = 25, mech_kwargs...)

data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
global const ROTATE = Ref{Float64}(0.0)
# @benchmark loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
# @profiler loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[loss(:box, pairs0, data0 + [i;zeros(55)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[loss(:box, pairs0, data0 + [0;i;zeros(54)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: L-BFGS:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################
using Optim

function termination_callback(opt; value_tol=1e-4)
	return opt.value < value_tol
end

solver = LBFGS(;m=100,
        alphaguess=Optim.LineSearches.InitialStatic(),
        linesearch=Optim.LineSearches.HagerZhang(),
        P = 1e1*I(25),
		)

lower = [0.00,
	+0.05, +0.05, -1.00,
	+0.05, -1.00, -1.00,
	-1.00, +0.05, -1.00,
	-1.00, -1.00, -1.00,
	+0.05, +0.05, +0.05,
	+0.05, -1.00, +0.05,
	-1.00, +0.05, +0.05,
	-1.00, -1.00, +0.05]
upper = [0.80,
	+1.00, +1.00, -0.05,
	+1.00, -0.05, -0.05,
	-0.05, +1.00, -0.05,
	-0.05, -0.05, -0.05,
	+1.00, +1.00, +1.00,
	+1.00, -0.05, +1.00,
	-0.05, +1.00, +1.00,
	-0.05, -0.05, +1.00]

function d2data(d)
	cf = d[1]
	data = [cf; 0;0;0; +d[2:4];
			cf; 0;0;0; +d[5:7];
			cf; 0;0;0; +d[8:10];
			cf; 0;0;0; +d[11:13];
			cf; 0;0;0; +d[14:16];
			cf; 0;0;0; +d[17:19];
			cf; 0;0;0; +d[20:22];
			cf; 0;0;0; +d[23:25];
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(25))

function fg!(F,G,d)
	# do common computations here
	l, ∇ = loss(:box, pairs0, d2data(d), n_sample=175)
	if G != nothing
		G .= ∇d2data' * ∇
	end
	if F != nothing
		value = l
	    return value
	end
end

d0 = [0.40,
	+0.50, +0.50, -0.50,
	+0.50, -0.50, -0.50,
	-0.50, +0.50, -0.50,
	-0.50, -0.50, -0.50,
	+0.50, +0.50, +0.50,
	+0.50, -0.50, +0.50,
	-0.50, +0.50, +0.50,
	-0.50, -0.50, +0.50]

ROTATE[] = 0.0
optimize(Optim.only_fg!(fg!), lower, upper, d0, Fminbox(solver),
	Optim.Options(
		callback = termination_callback,
		# allow_f_increases = true,
		show_trace = true),
	; inplace = false,
	)



# function callback(optim_state)
# 	# @show fieldnames(typeof(optim_state.metadata))
# 	# Main.SEED += 1/40
# 	return false
# end

# We can learn the coefficient of friction and the side dimenson of the cube
# form 15*0.75 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS


################################################################################
# Visualization
################################################################################
