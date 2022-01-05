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

Δt = 1/148
gscaled = -9.81*20

mech = getmechanism(:box, Δt=Δt, g=gscaled, cf=0.2, radius=0.00, side=0.50, mode=:box);
initialize!(mech, :box, x=[0,-1,1.], v=[0,2,1.], ω=[2,5,10.])
storage = simulate!(mech, 5.0, record=true,
    opts=InteriorPointOptions(btol=1e-6, rtol=1e-6, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)

################################################################################
# Generate & Save Dataset
################################################################################
generate_hardware_dataset(N=100, sleep_ratio=0.02)

################################################################################
# Load Dataset
################################################################################
params0, trajs0, pairs0 = open_dataset(:hardwarebox; N=100)

pairs0
params0
data0 = params0[:data]

################################################################################
# Optimization Objective: Evaluation & Gradient
################################################################################
clean_loss(:box, pairs0, data0, n_sample=250,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)


sss = [clean_loss(:box, pairs0[1:10], data0,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=i, g=gscaled)[1] for i = 0.001:0.0002:0.015]
plot(0.001:0.0002:0.015, log.(10, sss))

[clean_loss(:box, pairs0, data0 + [i;zeros(55)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)
	for i in Vector(-0.10:0.01:0.1)]
[clean_loss(:box, pairs0, data0 + [0;i;zeros(54)], opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)
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
	data = [cf; 0;0;0; +d[2]; +d[2]; -d[2];
			cf; 0;0;0; +d[2]; -d[2]; -d[2];
			cf; 0;0;0; -d[2]; +d[2]; -d[2];
			cf; 0;0;0; -d[2]; -d[2]; -d[2];
			cf; 0;0;0; +d[2]; +d[2]; +d[2];
			cf; 0;0;0; +d[2]; -d[2]; +d[2];
			cf; 0;0;0; -d[2]; +d[2]; +d[2];
			cf; 0;0;0; -d[2]; -d[2]; +d[2];
			]
	return data
end
∇d2data = ForwardDiff.jacobian(d -> d2data(d), zeros(2))

d0 = [0.40, 0.50]
	# +0.50, +0.50, -0.50,
	# +0.50, -0.50, -0.50,
	# -0.50, +0.50, -0.50,
	# -0.50, -0.50, -0.50,
	# +0.50, +0.50, +0.50,
	# +0.50, -0.50, +0.50,
	# -0.50, +0.50, +0.50,
	# -0.50, -0.50, +0.50]
lower = [0.00, 0.05]
	# +0.05, +0.05, -1.50,
	# +0.05, -1.50, -1.50,
	# -1.50, +0.05, -1.50,
	# -1.50, -1.50, -1.50,
	# +0.05, +0.05, +0.05,
	# +0.05, -1.50, +0.05,
	# -1.50, +0.05, +0.05,
	# -1.50, -1.50, +0.05]
upper = [0.80, 1.50]
	# +1.50, +1.50, -0.05,
	# +1.50, -0.05, -0.05,
	# -0.05, +1.50, -0.05,
	# -0.05, -0.05, -0.05,
	# +1.50, +1.50, +1.50,
	# +1.50, -0.05, +1.50,
	# -0.05, +1.50, +1.50,
	# -0.05, -0.05, +1.50]

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), n_sample=200, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=200, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=25, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-6)


clean_loss(:box, pairs0, d2data(dsol), n_sample=200,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)



f0(d0)
f0(dsol)
d0
dsol
# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS


using Plots; pyplot()
x=range(0.00,stop=0.80,length=20)
y=range(0.10,stop=1.50,length=20)
f(x,y) = log(10, clean_loss(:box, pairs0, d2data([x,y]), n_sample=200, Δt=Δt, g=gscaled,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1])
plot(x,y,f,st=:surface,camera=(20,60))






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

function f0(d; rot=0)
	return clean_loss(:box, pairs0, d2data(d), n_sample=1000, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))[1]
end

function fgH0(d; rot=0)
	f, g, H = clean_loss(:box, pairs0, d2data(d), n_sample=1000, Δt=Δt, g=gscaled, rot=rot, opts=InteriorPointOptions(btol=3e-4, rtol=3e-4))
	return f, ∇d2data' * g, ∇d2data' * H * ∇d2data
end
dsol, Dsol = quasi_newton_solve(f0, fgH0, d0, iter=80, gtol=1e-8, ftol=1e-6,
	lower=lower, upper=upper, reg=1e-6)


clean_loss(:box, pairs0, d2data(dsol), n_sample=200,
	opts=InteriorPointOptions(btol=3e-4, rtol=3e-4), Δt=Δt, g=gscaled)



f0(d0)
f0(dsol)
d0
dsol
# We can learn the coefficient of friction and the side dimenson of the cube
# form 35*0.40 seconds of recording. We use the simulator to evaluate the loss
# and its gradients by differentiating through the simulator. With gradient
# information we can use L-BFGS






using Polyhedra
# v = vrep(collect(permutations([0, 1, 2, 3])))
p4 = polyhedron()

v0 = vrep([[0, 0, 0], [0, 1, 1], [1/2, 1/2, 1], [1/2, 0/2, 3]])
v0 = vrep([dsol[1 + (i-1)*3 .+ (1:3)] for i = 1:8])
v0 = vrep([d0[1 + (i-1)*3 .+ (1:3)] for i = 1:8])
p0 = polyhedron(v0)
Polyhedra.Mesh(p0)
setobject!(vis, Polyhedra.Mesh(p0))


# Project that polyhedron down to 3 dimensions for visualization
v1 = [1, -1,  0,  0]
v2 = [1,  1, -2,  0]
v3 = [1,  1,  1, -3]
p3 = project(p4, [v1 v2 v3])

# Show the result
setobject!(vis, Polyhedra.Mesh(p3))












################################################################################
# Visualization
################################################################################

# Visualize Dataset with solution
datasol = d2data(dsol)
# datasol = d2data([0.1, 1.0])
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
	# datasol = d2data(d0)
	# datasol = d2data(lower)
	# datasol = d2data(upper)
	# datasol = data0
	for i = 1:8
		datasol[(i-1)*7 + 4] = 0.05
	end
	set_simulator_data!(mech, datasol)
	visualize(mech, trajs0[i], vis=vis, show_contact = true)
	sleep(1.5)
end
dsol
Dsol

mech.Δt


d2data(lower)
d2data(upper)




z1 = pairs0[1][1]
z2 = pairs0[2][1]
z3 = pairs0[3][1]
z4 = pairs0[4][1]
z5 = pairs0[5][1]
z6 = pairs0[6][1]
z7 = pairs0[7][1]

x1, v15, q1, ϕ15 = unpackMaxState(z1, 1)
x2, v25, q2, ϕ25 = unpackMaxState(z2, 1)
x3, v35, q3, ϕ35 = unpackMaxState(z3, 1)
x4, v45, q4, ϕ45 = unpackMaxState(z4, 1)
x5, v55, q5, ϕ55 = unpackMaxState(z5, 1)
x6, v65, q6, ϕ65 = unpackMaxState(z6, 1)
x7, v75, q7, ϕ75 = unpackMaxState(z7, 1)


v25b = (x2 - x1) / Δt
v35b = (x3 - x2) / Δt
v45b = (x4 - x3) / Δt
v55b = (x5 - x4) / Δt
v65b = (x6 - x5) / Δt
v75b = (x7 - x6) / Δt

(v35b - v25b) / (Δt * 9.81)
(v45b - v35b) / (Δt * 9.81)
(v55b - v45b) / (Δt * 9.81)
(v65b - v55b) / (Δt * 9.81)
(v75b - v65b) / (Δt * 9.81)

9.81*20

(v25 - v15) / Δt
(v35 - v25) / Δt
(v45 - v35) / Δt
(v55 - v45) / Δt
(v65 - v55) / Δt
(v75 - v65) / Δt

(x2 - (x1 + Δt * v15)) ./ (Δt * v15)
(x2 - (x1 + Δt * v25)) ./ (Δt * v25)
(x2 - (x1 + Δt * (v15+v25)/2)) ./ (Δt * (v15+v25)/2)

(x3 - (x2 + Δt * v25)) ./ (Δt * v25)
(x3 - (x2 + Δt * v35)) ./ (Δt * v35)
(x3 - (x2 + Δt * (v25+v35)/2)) ./ (Δt * (v25+v35)/2)

(x4 - (x3 + Δt * v35)) ./ (Δt * v35)
(x4 - (x3 + Δt * v45)) ./ (Δt * v45)
(x4 - (x3 + Δt * (v35+v45)/2)) ./ (Δt * (v35+v45)/2)

(x5 - (x4 + Δt * v45)) ./ (Δt * v45)
(x5 - (x4 + Δt * v55)) ./ (Δt * v55)
(x5 - (x4 + Δt * (v45+v55)/2)) ./ (Δt * (v45+v55)/2)

plot(hcat(x1, x2, x3, x4, x5)')
plot(hcat([p[1][1:3] for p in pairs0[1:15]]...)' )


mean([(p[2][4:6] - p[1][4:6])[3]/ (Δt * -9.81) for p in pairs0[1:15]])
median([(p[2][4:6] - p[1][4:6])[3]/ (Δt * -9.81) for p in pairs0[1:15]])
plot([(p[2][4:6] - p[1][4:6])[3]/ (Δt * -9.81) for p in pairs0[1:300]])

Δt = 1/148
-9.81 * Δt
(v25 - v15)/ Δt
(v35 - v25)/ Δt
(v45 - v35)/ Δt
(v55 - v45)/ Δt

v15



v25



v35



v45



v55

Δt = 1/148
function ll(Δt)
	x2t = x1 + Δt*v15
	v25t = v15 + Δt*[0,0,18.7 * -9.81]
	ϕ25t = ϕ15 + Δt*[0,0,-9.81]
	q2t = getq3(q1, SVector{3}(ϕ15), Δt)
	Δx2 = norm(x2t - x2)
	Δv25 = norm(v25t - v25)
	Δq2 = norm(vector(q2t) - vector(q2))
	return [Δx2, Δv25, Δq2]
end

plt = plot()
Δx2 = hcat([log.(10, ll(t)) for t in 0.0001:0.0001:0.01]...)[1,:]
Δv25 = hcat([log.(10, ll(t)) for t in 0.0001:0.0001:0.01]...)[2,:]
Δq2 = hcat([log.(10, ll(t)) for t in 0.0001:0.0001:0.01]...)[3,:]
plot!(plt, 0.0001:0.0001:0.01, Δx2, label="Δx2")
plot!(plt, 0.0001:0.0001:0.01, Δv25, label="Δv25")
plot!(plt, 0.0001:0.0001:0.01, Δq2, label="Δq2")
scatter!(plt, [1/148], [0])
x2 - (x1 + Δt * ((v25 + v15) / 2))



a = 10
a = 10
a = 10
a = 10
a = 10
a = 10
