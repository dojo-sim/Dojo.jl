using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# ## Setup
using Dojo 
include("methods/simulation.jl")
include("methods/process.jl")
include("methods/benchmark.jl")

# ## Visualizer
vis = Visualizer()
open(vis)

################################################################################
## Astronaut
################################################################################
## humanoid
## multiple bodies
## no initial linear velocities
## no gravity
## no spring and damper
## random control
################################################################################

function ctrl!(mechanism, t)
	nu = input_dimension(mechanism)
	u = 0.3 * [zeros(6); mechanism.timestep; zeros(nu-7)]
	set_input!(mechanism, u)
	return
end

# ## Initialize
Random.seed!(0)
ϵ = 1.0e-14

# ## Mechanism
mech = get_mechanism(:humanoid, 
	timestep=0.01, 
	gravity=0.0, 
	spring=0.0, 
	damper=0.0, 
	contact=false)

# ## Simulate
initialize!(mech, :humanoid)
storage = simulate!(mech, 1.0, ctrl!, 
	record=true, 
	opts=SolverOptions(rtol=ϵ, btol=ϵ))

# ## Visualize
visualize(mech, storage, 
	vis=vis)

################################################################################
## Momentum
################################################################################

# ## Simulate
timestep= [0.10,]
mech, storage, tcompute = astronaut_simulation(;
	Nsim=1, 
	timestep=1e-1, 
	tsim=10.0, 
	tctrl=1.0, 
	seed=0, 
	control_amplitude=0.05)

# ## Visualize
visualize(mech, storage[1], vis=vis)

# ## Benchmark
timestep= [0.10, 0.03, 0.01, 0.003, 0.001]
speed_dj, ener_dj = benchmark_energy(timestep; Nsim=1, tsim=2.0, tctrl=1.0, seed=0, ϵ=1e-14)
speed_dj, plin_dj, pang_dj = benchmark_momentum(timestep; Nsim=1, tsim=1.0, seed=0, ϵ=1e-14)

# ## Saving results
path = joinpath(@__DIR__, "results/astronaut_momentum.jld2")
jldsave(path, speed=speed_dj, plin=plin_dj, pang=pang_dj)

# ## Recover results
file = jldopen(joinpath(@__DIR__, "../MuJoCoBenchmark.jl/results/astronaut_momentum.jld2"))
speed_mj, plin_mj, pang_mj, ener_mj = file["speed"], file["plin"], file["pang"], file["ener"]
file = jldopen(joinpath(@__DIR__, "results/astronaut_momentum.jld2"))
speed_dj, plin_dj, pang_dj = file["speed"], file["plin"], file["pang"]

plt = plot(layout=(1,2), legend=false, xlims=(1e-2,2e4), ylims=(1e-15,1e2))
plot!(plt[1,1], speed_dj, plin_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed_dj, pang_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed_dj, plin_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed_dj, pang_dj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)

plot!(plt[1,1], speed_mj, plin_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed_mj, pang_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed_mj, plin_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed_mj, pang_mj, axis=:log, yaxis=:flip, linewidth=3, markersize=6)

################################################################################
## Energy
################################################################################

# ## Simulate
timestep= [0.10, 0.02, 0.005]
ener_traj = []
for i = 1:3
	mech, storage, tcompute = astronaut_simulation(;
		Nsim=1, 
		timestep=timestep[i], 
		tsim=100.0, 
		tctrl=1.0, 
		seed=0, 
		control_amplitude=0.02)
	stride = Int(floor(1 / timestep[i]))
	ener_t = [kinetic_energy(mech, storage[1], i) for i in stride+2:stride:length(storage[1])]
	push!(ener_traj, ener_t .- ener_t[1])
	@show i
end

# ## Saving results
path = joinpath(@__DIR__, "results/astronaut_energy.jld2")
jldsave(path, ener_traj=ener_traj)

# ## Recover results
file = jldopen(joinpath(@__DIR__, "../MuJoCoBenchmark.jl/results/astronaut_energy.jld2"))
ener_traj_mj = file["ener_traj"]
file = jldopen(joinpath(@__DIR__, "results/astronaut_energy.jld2"))
ener_traj_dj = file["ener_traj"]

# ## Plot
plt = plot(layout=(3,1), legend=false, xlims=(0,100))
for i = 1:3
	plot!(plt[i,1],
		range(0,stop=100,length=length(ener_traj_mj[i])),
		ener_traj_mj[i], linewidth=3, markersize=6)
	plot!(plt[i,1],
		range(0,stop=100,length=length(ener_traj_dj[i])),
		ener_traj_dj[i], linewidth=3, markersize=6)
end
display(plt)



