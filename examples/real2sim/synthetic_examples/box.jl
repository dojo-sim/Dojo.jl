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

mech = getmechanism(:box, Δt=0.05, g=-9.81, cf=0.2, radius=0.00, side=0.50, mode=:box);
initialize!(mech, :box, x=[0,-1,1.], v=[0,2,1.], ω=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=false)


################################################################################
# Generate & Save Dataset
################################################################################
init_kwargs = Dict(:xlims => [[0,0,0.2], [1,1,0.4]],
				   :vlims => [-2ones(3), [2,2,-1.]],
				   :ωlims => [-6ones(3), 6ones(3)])
mech_kwargs = Dict(:cf => 0.1, :radius => 0.0, :side => 0.5)
generate_dataset(:box, H=0.40, N=35,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4),
	init_kwargs=init_kwargs,
	mech_kwargs=mech_kwargs,
	sleep_ratio=0.01,
	show_contact=false)


################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:box; N = 35, mech_kwargs...)

data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
# @benchmark clean_loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
# @profiler clean_loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
clean_loss(:box, pairs0, data0, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))

[clean_loss(:box, pairs0, data0 + [i;zeros(55)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:box, pairs0, data0 + [0;i;zeros(54)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: Quasi Newton:
# We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################


include("../quasi_newton.jl")
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

d0 = [0.40,
	+0.50, +0.50, -0.50,
	+0.50, -0.50, -0.50,
	-0.50, +0.50, -0.50,
	-0.50, -0.50, -0.50,
	+0.50, +0.50, +0.50,
	+0.50, -0.50, +0.50,
	-0.50, +0.50, +0.50,
	-0.50, -0.50, +0.50]
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

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), n_sample=35, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=35, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol = quasi_newton_solve(f0, fgH0, d0, iter=50, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

d0
f0(d0)
f0(dsol)
norm(data0 - d2data(dsol), Inf)

# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS


################################################################################
# Visualization
################################################################################

include("../quasi_newton.jl")
function d2data(d)
	data = [d[1]; 0;0;0; +d[2:4];
			d[5]; 0;0;0; +d[6:8];
			d[9]; 0;0;0; +d[10:12];
			d[13]; 0;0;0; +d[14:16];
			d[17]; 0;0;0; +d[18:20];
			d[21]; 0;0;0; +d[22:24];
			d[25]; 0;0;0; +d[26:28];
			d[29]; 0;0;0; +d[30:32];
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(32))

d0 = [
	0.40, +0.50, +0.50, -0.50,
	0.40, +0.50, -0.50, -0.50,
	0.40, -0.50, +0.50, -0.50,
	0.40, -0.50, -0.50, -0.50,
	0.40, +0.50, +0.50, +0.50,
	0.40, +0.50, -0.50, +0.50,
	0.40, -0.50, +0.50, +0.50,
	0.40, -0.50, -0.50, +0.50]
lower = [
	0.00, +0.05, +0.05, -1.00,
	0.00, +0.05, -1.00, -1.00,
	0.00, -1.00, +0.05, -1.00,
	0.00, -1.00, -1.00, -1.00,
	0.00, +0.05, +0.05, +0.05,
	0.00, +0.05, -1.00, +0.05,
	0.00, -1.00, +0.05, +0.05,
	0.00, -1.00, -1.00, +0.05]
upper = [
	0.80, +1.00, +1.00, -0.05,
	0.80, +1.00, -0.05, -0.05,
	0.80, -0.05, +1.00, -0.05,
	0.80, -0.05, -0.05, -0.05,
	0.80, +1.00, +1.00, +1.00,
	0.80, +1.00, -0.05, +1.00,
	0.80, -0.05, +1.00, +1.00,
	0.80, -0.05, -0.05, +1.00]

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), n_sample=35, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=35, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end

dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=50, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-9)

d0
f0(d0)
f0(dsol)
norm(data0 - d2data(dsol), Inf)

# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS


################################################################################
# Visualization
################################################################################

# Visualize Dataset with solution
datasol = ∇d2data * dsol
for i = 1:8
	datasol[(i-1)*7 + 4] = 0.05
end
set_simulator_data!(mech, datasol)
for traj in trajs0
	visualize(mech, traj, vis=vis, show_contact = true)
	sleep(0.5)
end

# Visualize Progress made during 'learning'
for (i,dsol) in enumerate(Dsol[1:6])
	datasol = ∇d2data * dsol
	for i = 1:8
		datasol[(i-1)*7 + 4] = 0.05
	end
	set_simulator_data!(mech, datasol)
	visualize(mech, trajs0[i], vis=vis, show_contact = true)
	sleep(1.5)
end
