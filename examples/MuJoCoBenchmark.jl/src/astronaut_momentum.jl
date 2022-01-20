# Utils
function module_dir()
    return joinpath(@__DIR__, "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using LinearAlgebra
using LyceumBase
using LyceumMuJoCo
using LyceumMuJoCoViz
using MuJoCo
using Plots
using JLD2
using Random

include("methods.jl")

function mj_astronaut_simulation(; Δt=0.01, tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.05, record_traj::Bool=true, record_ener::Bool=false)

	jm, jd, mjsim = mj_model("astronaut.xml", Δt=Δt)
    N = Int(floor(tsim/Δt))
	tcompute = 0.0
	ener0 = 0.0
	ener_traj = []

	traj = zeros(83,0)
	for i = 1:N
		jd.ctrl .= (i*Δt < tctrl) * control_amplitude * rand(21)
		tcompute += @elapsed mj_step(jm, jd)
		record_traj && (traj = hcat(traj, deepcopy(getstate(mjsim))))
		if !(i*Δt < tctrl) && ((i-1)*Δt < tctrl)
			# first step after the end of control inputs
			ener0 = energy(jm, jd)
		end
		if record_ener && !(i*Δt < tctrl) && (i % Int(floor(1/Δt)) == 0)
			push!(ener_traj, energy(jm, jd) - ener0)
		end
	end
	fail = norm(jd.qvel,Inf) < 1e-10 # when Mujoco fails it returns a zero qvel
	ener1 = energy(jm, jd)
	ener = ener1 - ener0
	plin = norm(linear_momentum(jm, jd))
	pang = norm(angular_momentum(jm, jd))
    return tcompute, plin, pang, ener, fail, traj, ener_traj
end

function mj_batch_astronaut_simulation(; Nsim=1, Δt=0.01, tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.05, record_traj::Bool=false)
	speed = 0.0
	plin = 0.0
	pang = 0.0
	ener = 0.0
	fail = 0.0
	for i = 1:Nsim
		tcompute, pl, pa, en, fa, _ = mj_astronaut_simulation(; Δt=Δt, tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ,
			control_amplitude=control_amplitude, record_traj=record_traj)
		speed += tsim / tcompute
		plin += pl
		pang += pa
		ener += abs(en)
		fail += fa
	end
	speed /= Nsim
	plin /= Nsim
	pang /= Nsim
	ener /= Nsim
	fail /= Nsim
    return speed, plin, pang, ener, fail
end

function mj_benchmark_momentum(Δt; Nsim=1, tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.05, record_traj::Bool=false)
	N = length(Δt)
	speed = zeros(N)
	plin = zeros(N)
	pang = zeros(N)
	ener = zeros(N)
	fail = zeros(N)
	for i = 1:N
		@show Δt[i]
		speed[i], plin[i], pang[i], ener[i], fail[i] = mj_batch_astronaut_simulation(; Nsim=Nsim, Δt=Δt[i],
			tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ, control_amplitude=control_amplitude, record_traj=record_traj)
	end
    return speed, plin, pang, ener, fail
end


speed, plin, pang, Δener, fail, traj, ener_traj = mj_astronaut_simulation(tsim=1, tctrl=1.0,
	Δt=0.1, control_amplitude=0.02)
speed, plin, pang, Δener, fail, traj, ener_traj = mj_astronaut_simulation(tsim=100, tctrl=1.0,
	Δt=0.1, control_amplitude=0.02)
LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

################################################################################
# momentum
################################################################################
# mj_batch_astronaut_simulation(Nsim=10)
Δt = [0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003]
speed, plin, pang, ener, fail = mj_benchmark_momentum(Δt, Nsim=3, control_amplitude=0.05)

# Saving results
# path = joinpath(module_dir(), "results/astronaut_momentum.jld2")
# jldsave(path, speed=speed, plin=plin, pang=pang)

plt = plot(layout=(1,2), legend=false, xlims=(1e-2,1e6), ylims=(1e-15,1e6))
plot!(plt[1,1], speed, plin, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed, pang, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed, plin, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed, pang, axis=:log, yaxis=:flip, linewidth=3, markersize=6)

################################################################################
# energy
################################################################################
Δt = [0.1, 0.01, 0.001]
ener_traj = []
for i = 1:3
	speed, plin, pang, Δener, fail, traj, ener_t = mj_astronaut_simulation(tsim=100, tctrl=1.0,
		Δt=Δt[i], control_amplitude=0.02, record_traj=false, record_ener=true)
	push!(ener_traj, ener_t)
end

plt = plot(layout=(3,1), legend=false)#, xlims=(0,100))
for i = 1:3
	plot!(plt[i,1],
		range(0,stop=100,length=length(ener_traj[i])),
		ener_traj[i], linewidth=3, markersize=6)
end
display(plt)


# Saving results
path = joinpath(module_dir(), "results/astronaut_energy.jld2")
jldsave(path, ener_traj=ener_traj)




################################################################################
# Demo
################################################################################

Δt = 0.01
jm, jd, mjsim = mj_model("astronaut.xml", Δt=0.01)
tsim = 1.0
N = Int(floor(tsim/Δt))

# jd.qpos[3] += 10.0
# Random.seed!(100)
# jd.qvel .= 1.0 * (rand(27) .- 0.5)

traj = zeros(83, 0)
E = zeros(N)
Plin = zeros(N)
Pang = zeros(N)
tcompute = 0.0
for i = 1:N
    jd.ctrl .= 0.3 * [0; 1; zeros(19)]
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
    E[i] = energy(jm, jd)
	Plin[i] = norm(linear_momentum(jm, jd))
	# Pang[i] = angular_momentum(jm, jd)
end
tsim / tcompute
plot(E)
plot(Plin)
plot(Pang)

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel
