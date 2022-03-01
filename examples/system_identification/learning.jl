using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# ## Setup
using Dojo
using Random
using LinearAlgebra 
using Plots

include("methods/utils.jl")
include("methods/data_processing.jl")
include("methods/quasi_newton.jl")
include("hardware_examples/methods.jl")
include("hardware_examples/visualization.jl")

# ## Open visualizer
vis= Visualizer()
open(vis)

# ## Mechanism
S = 7
timestep= 1/148 * S
gscaled = -9.81*20
mech = get_mechanism(:box, 
	timestep=timestep, 
	gravity=gravityscaled, 
	friction_coefficient=0.2, 
	radius=0.00, 
	side=0.50, 
	mode=:box);

# ## Simualate
initialize!(mech, :box, 
	x=[0.0, -1.0, 1.0], 
	v=[0.0, 2.0, 1.0], 
	ω=[2.0, 5.0, 10.0])
storage = simulate!(mech, 5.0, 
	record=true,
    opts=SolverOptions(btol=1e-6, rtol=1e-6, verbose=false))

# ## Visualize
visualize(mech, storage, vis=vis, show_contact=true)

################################################################################
## Generate & Save Dataset
################################################################################
generate_hardware_dataset(N=400, sleep_ratio=0.0002, S=S)
generate_hardware_dataset(N=400, sleep_ratio=0.0002, S=1)

################################################################################
## Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:hardwarebox; N=400, S=S)
params1, trajs1, pairs1 = open_dataset(:hardwarebox; N=400, S=1)

pairs0
params0
data0 = params0[:data]

################################################################################
## Optimization Objective: Evaluation & Gradient
################################################################################
clean_loss(:box, pairs0, data0, 
	n_sample=250,
	opts=SolverOptions(btol=3e-4, rtol=3e-4), 
	timestep=timestep, 
	gravity=gravityscaled)

[clean_loss(:box, pairs0, data0 + [i;zeros(55)],
	opts=SolverOptions(btol=3e-4, rtol=3e-4), 
	timestep=timestep, 
	gravity=gravityscaled)
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:box, pairs0, data0 + [0;i;zeros(54)],
	opts=SolverOptions(btol=3e-4, rtol=3e-4), 
	timestep=timestep, 
	gravity=gravityscaled)
	for i in Vector(-0.10:0.01:0.1)]

plot(hcat([p[1][1:3] for p in pairs0]...)')
plot(hcat([p[1][4:6] for p in pairs0]...)')
plot(hcat([p[1][7:10] for p in pairs0]...)')
plot(hcat([p[1][11:13] for p in pairs0]...)')

################################################################################
## Optimization Algorithm: Quasi Newton:
## We learn a single coefficient of friction and a side length = 2 params
################################################################################
function d2data(d)
	friction_coefficient = d[1]
	side = d[2]
	data = [friction_coefficient; 0.0; 0.0; 0.0; +side; +side; -side;
			friction_coefficient; 0.0; 0.0; 0.0; +side; -side; -side;
			friction_coefficient; 0.0; 0.0; 0.0; -side; +side; -side;
			friction_coefficient; 0.0; 0.0; 0.0; -side; -side; -side;
			friction_coefficient; 0.0; 0.0; 0.0; +side; +side; +side;
			friction_coefficient; 0.0; 0.0; 0.0; +side; -side; +side;
			friction_coefficient; 0.0; 0.0; 0.0; -side; +side; +side;
			friction_coefficient; 0.0; 0.0; 0.0; -side; -side; +side;
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 0.50]
lower = [0.00, 0.05]
upper = [0.80, 1.50]

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), 
		n_sample=200, 
		timestep=timestep, 
		gravity=gravityscaled, 
		rot=rot, 
		opts=SolverOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; 
	rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), 
		n_sample=200, 
		timestep=timestep, 
		gravity=gravityscaled, 
		rot=rot, 
		opts=SolverOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, 
	iter=100, 
	gtol=1e-8, 
	ftol=1e-6,
	lower=lower, 
	upper=upper, 
	reg=1e-6, 
	Δrot=20)

f0(d0)
f0(dsol)
d0
dsol #[0.171, 0.933]
## We can learn the coefficient of friction and the side dimenson of the cube
## form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
## and its gradients by differentiating through the simulator. With gradient
## information we can use L-BFGS

################################################################################
## Visualization
################################################################################

using Plots; pyplot()
x=range(0.00,stop=0.80,length=10)
y=range(0.10,stop=1.50,length=10)
f(x,y) = log(10, clean_loss(:box, pairs0, d2data([x,y]), n_sample=50, timestep=timestep, gravity=gravityscaled,
	opts=SolverOptions(btol=3e-4, rtol=3e-4))[1])
plot(x,y,f,st=:surface,camera=(25,75))


################################################################################
## Optimization Algorithm: Quasi Newton:
## We learn a single coefficient of friction and a 8 contact locations [x,y,z] -> 25 params in total
################################################################################
function d2data(d)
	friction_coefficient = d[1]
	data = [friction_coefficient; 0;0;0; +d[2:4];
			friction_coefficient; 0;0;0; +d[5:7];
			friction_coefficient; 0;0;0; +d[8:10];
			friction_coefficient; 0;0;0; +d[11:13];
			friction_coefficient; 0;0;0; +d[14:16];
			friction_coefficient; 0;0;0; +d[17:19];
			friction_coefficient; 0;0;0; +d[20:22];
			friction_coefficient; 0;0;0; +d[23:25];
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
## d0 = [0.40,
## 	+0.50, +0.50, -0.50,
## 	+0.50, -0.50, -0.50,
## 	-0.50, +0.50, -0.50,
## 	-0.50, -0.50, -0.50,
## 	+0.50, +0.50, +0.50,
## 	+0.50, -0.50, +0.50,
## 	-0.50, +0.50, +0.50,
## 	-0.50, -0.50, +0.50]
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
	return clean_loss(:box, pairs0, d2data(d), n_sample=n_sample, timestep=timestep, gravity=gravityscaled, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0, n_sample=50)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=n_sample, timestep=timestep, gravity=gravityscaled, rot=rot, opts=SolverOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=50, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e+1, Δrot=500, n_sample0=20, Δn_sample=5, n_sample_max=1000)

f0(d0)
f0(dsol)
d0
dsol
## We can learn the coefficient of friction and the side dimenson of the cube
## form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
## and its gradients by differentiating through the simulator. With gradient
## information we can use L-BFGS

jldsave(joinpath("..", "results", "sol.jld2"), dsol=dsol, Dsol=Dsol)


################################################################################
## Visualization
################################################################################

# ## Visualize Dataset with solution
datasol = d2data(dsol)
for i = 1:8
	datasol[(i-1)*7 + 4] = 0.05
end
set_simulator_data!(mech, datasol)
for traj in trajs0[1:10]
	visualize(mech, traj, vis=vis, show_contact = true)
	sleep(0.5)
end

# ## Visualize Progress made during 'learning'
for (i,dsol) in enumerate(Dsol[1:20])
	datasol = d2data(dsol)
	for i = 1:8
		datasol[(i-1)*7 + 4] = 0.05
	end
	set_simulator_data!(mech, datasol)
	visualize(mech, trajs0[i], vis=vis, show_contact = true)
	sleep(1.5)
end
