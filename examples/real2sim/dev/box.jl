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
visualize(mech, storage, vis=vis, show_contact = true)


dbest = [0.09902443899194172,
	0.24983056763668757, 0.2501325129483943, -0.25042971065786357,
	0.25000879472930276, -0.25051893328360614, -0.2501886552502975,
	-0.24996123416616736, 0.25013736791624774, -0.25006562608076943,
	-0.24431161099368243, -0.2500961647900991, -0.2502248998529076,
	0.2500992241490191, 0.25005706863989413, 0.25006535845090355,
	0.2497660464529427, -0.2503818791472568, 0.25004117849998964,
	-0.24991169274364852, 0.2500621347605857, 0.24994574048892784,
	-0.25023908162014147, -0.2499294194836169, 0.2499100511808172]

dsol = [0.08411049246523367, 0.1004263971744767, 0.6949490819623005,
	-0.22376939359194972, 0.2577953843534969, -0.249872662027289,
	-0.2504716235835994, -0.23307586778080344, 0.6354201175987965,
	-0.24942840292413965, -0.25076835774075035, -0.24992475320312135,
	-0.25006651351196774, 0.257147320863724, 0.7743336515556084,
	0.24746607793048037, 0.24674202858214214, -0.2500290803343864,
	0.2502303239464445, -0.2680786838158817, 0.5830544824156076,
	0.2391490420345988, -0.24503110167092726, -0.2501244181194456, 0.2505384341205768]

datasol = ∇d2data * dsol
datasol = ∇d2data * dbest
for i = 1:8
	datasol[(i-1)*7 + 4] = 0.05
end
set_simulator_data!(mech, datasol)
for traj in trajs0
	visualize(mech, traj, vis=vis, show_contact = true)
	sleep(0.5)
end


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

solver = LBFGS(;m=10,
# solver = BFGS(
        alphaguess=Optim.LineSearches.InitialStatic(),
        linesearch=Optim.LineSearches.HagerZhang(),
        P = 1e2*I(25),
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
	l, ∇ = loss(:box, pairs0, d2data(d), n_sample=50)
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
