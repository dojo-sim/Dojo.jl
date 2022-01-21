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

S = 7
Δt = 1/148 * S
gscaled = -9.81*20

mech = getmechanism(:box, Δt=Δt, g=gscaled, cf=0.2, radius=0.00, side=0.50, mode=:box);
initialize!(mech, :box, x=[0,-1,1.], v=[0,2,1.], ω=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)

################################################################################
# Generate & Save Dataset
################################################################################
generate_hardware_dataset(N=400, sleep_ratio=0.0002, S=S)
generate_hardware_dataset(N=400, sleep_ratio=0.0002, S=1)

################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:hardwarebox; N=400, S=S)
params1, trajs1, pairs1 = open_dataset(:hardwarebox; N=400, S=1)

pairs0
params0
data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
clean_loss(:box, pairs0, data0, n_sample=250,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)

[clean_loss(:box, pairs0, data0 + [i;zeros(55)],
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:box, pairs0, data0 + [0;i;zeros(54)],
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')


################################################################################
# Optimization Algorithm: Quasi Newton:
# We learn a single coefficient of friction and a side length = 2 params
################################################################################
include("../quasi_newton.jl")
function d2data(d)
	cf = d[1]
	side = d[2]
	data = [cf; 0;0;0; +side; +side; -side;
			cf; 0;0;0; +side; -side; -side;
			cf; 0;0;0; -side; +side; -side;
			cf; 0;0;0; -side; -side; -side;
			cf; 0;0;0; +side; +side; +side;
			cf; 0;0;0; +side; -side; +side;
			cf; 0;0;0; -side; +side; +side;
			cf; 0;0;0; -side; -side; +side;
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 0.50]
lower = [0.00, 0.05]
upper = [0.80, 1.50]

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), n_sample=200, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=200, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=100, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-6, Δrot=20)

f0(d0)
f0(dsol)
d0
dsol #[0.171, 0.933]
# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

################################################################################
# Visualization
################################################################################

using Plots; pyplot()
x=range(0.00,stop=0.80,length=10)
y=range(0.10,stop=1.50,length=10)
f(x,y) = log(10, clean_loss(:box, pairs0, d2data([x,y]), n_sample=50, Δt=Δt, g=gscaled,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1])
plot(x,y,f,st=:surface,camera=(25,75))


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
	+1.30, +1.40, -1.50,
	+1.50, -1.30, -1.10,
	-1.50, +1.30, -1.10,
	-1.50, -1.40, -1.30,
	+1.10, +1.20, +1.50,
	+1.20, -1.40, +1.50,
	-1.50, +1.40, +1.30,
	-1.10, -1.40, +1.50]
# d0 = [0.40,
# 	+0.50, +0.50, -0.50,
# 	+0.50, -0.50, -0.50,
# 	-0.50, +0.50, -0.50,
# 	-0.50, -0.50, -0.50,
# 	+0.50, +0.50, +0.50,
# 	+0.50, -0.50, +0.50,
# 	-0.50, +0.50, +0.50,
# 	-0.50, -0.50, +0.50]
lower = [0.00,
	+0.05, +0.05, -1.50,
	+0.05, -1.50, -1.50,
	-1.50, +0.05, -1.50,
	-1.50, -1.50, -1.50,
	+0.05, +0.05, +0.05,
	+0.05, -1.50, +0.05,
	-1.50, +0.05, +0.05,
	-1.50, -1.50, +0.05]
upper = [0.80,
	+1.50, +1.50, -0.05,
	+1.50, -0.05, -0.05,
	-0.05, +1.50, -0.05,
	-0.05, -0.05, -0.05,
	+1.50, +1.50, +1.50,
	+1.50, -0.05, +1.50,
	-0.05, +1.50, +1.50,
	-0.05, -0.05, +1.50]

function f0(d; rot=0, n_sample=50)
	return clean_loss(:box, pairs0, d2data(d), n_sample=n_sample, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0, n_sample=50)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=n_sample, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=50, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e+1, Δrot=500, n_sample0=20, Δn_sample=5, n_sample_max=1000)

f0(d0)
f0(dsol)
d0
dsol
# x:[" 3.5e-1", " 1.0e+0", " 8.2e-1", "-1.0e+0", " 8.1e-1", "-1.0e+0", "-10.0e-1", "-8.0e-1", " 9.7e-1", "-1.0e+0", "-9.4e-1", "-1.1e+0", "-9.4e-1", " 9.4e-1", " 1.0e+0", " 1.0e+0", " 9.3e-1", "-1.0e+0", " 1.0e+0", "-9.4e-1", " 1.1e+0", " 9.4e-1", "-9.2e-1", "-1.1e+0", " 8.9e-1"]
# x:[" 3.2e-1", " 7.6e-1", " 5.9e-1", "-3.2e-1", " 7.6e-1", "-5.8e-1", "-5.9e-1", "-6.6e-1", " 4.7e-1", "-2.5e-1", "-7.5e-1", "-7.1e-1", "-4.6e-1", " 4.8e-1", " 6.7e-1", " 6.1e-1", " 6.0e-1", "-6.8e-1", " 2.7e-1", "-6.2e-1", " 6.1e-1", " 5.1e-1", "-3.5e-1", "-4.2e-1", " 6.1e-1"]

# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS

jldsave(joinpath(module_dir(), "examples", "real2sim", "hardware_examples", "sol.jld2"), dsol=dsol, Dsol=Dsol)


################################################################################
# Visualization
################################################################################

# Visualize Dataset with solution
datasol = d2data(dsol)
for i = 1:8
	datasol[(i-1)*7 + 4] = 0.05
end
set_simulator_data!(mech, datasol)
for traj in trajs0[1:10]
	visualize(mech, traj, vis=vis, show_contact = true)
	sleep(0.5)
end

# Visualize Progress made during 'learning'
for (i,dsol) in enumerate(Dsol[1:20])
	datasol = d2data(dsol)
	for i = 1:8
		datasol[(i-1)*7 + 4] = 0.05
	end
	set_simulator_data!(mech, datasol)
	visualize(mech, trajs0[i], vis=vis, show_contact = true)
	sleep(1.5)
end






# # Open visualizer
vis = Visualizer()
open(vis)

mech = getmechanism(:box, Δt=Δt/S, g=gscaled, cf=Dsol[end][1], radius=0.00, side=2.0, mode=:box);
set_simulator_data!(mech, d2data(Dsol[end]))
id = 7#4,6,7,8
traj_truth = trajs1[id]
x2 = traj_truth.x[1][1] - [0,0,2.0]/2
v15 = traj_truth.v[1][1]
q2 = traj_truth.q[1][1]
ϕ15 = traj_truth.ω[1][1]

initialize!(mech, :box, x=x2, v=v15, q=q2, ω=ϕ15)
traj_sim = simulate!(mech, 0.80, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))

cube_sim_v_truth(Dsol[end], traj_truth, traj_sim, vis=vis,
	transparency_truth=1.0,
	fps=Int(floor(1/mech.Δt)), b0=0.0, b1=0.0)

cube_ghost_sim_v_truth(Dsol[end], traj_truth, traj_sim, vis=vis, transparency_truth=1.0)
