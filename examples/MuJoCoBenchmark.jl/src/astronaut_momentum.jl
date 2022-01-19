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
		control_amplitude=0.1)

	jm, jd, mjsim = mj_model(:astronaut, Δt=Δt)
    N = Int(floor(tsim/Δt))
	plin = 0.0
	pang = 0.0
	ener = 0.0
	tcompute = 0.0

	traj = zeros(83,0)
	for i = 1:N
		jd.ctrl .= (i*Δt < tctrl) * control_amplitude * Δt * rand(21)
		tcompute += @elapsed mj_step(jm, jd)
		traj = hcat(traj, deepcopy(getstate(mjsim)))
	end
	ener = energy(jm, jd)
	plin = norm(linear_momentum(jm, jd))
	# pang = norm(angular_momentum(jm, jd))
    return tcompute, plin, pang, ener
end

function mj_batch_astronaut_simulation(; Nsim=1, Δt=0.01, tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.1)
	speed = 0.0
	plin = 0.0
	pang = 0.0
	ener = 0.0
	for i = 1:Nsim
		tcompute, pl, pa, en = mj_astronaut_simulation(; Δt=Δt, tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ,
			control_amplitude=control_amplitude)
		speed += tsim / tcompute
		plin += pl
		pang += pa
		ener += en
	end
	speed /= Nsim
	plin /= Nsim
	pang /= Nsim
	ener /= Nsim
    return speed, plin, pang, ener
end

function mj_benchmark_momentum(Δt; Nsim=1, tsim=1.0, tctrl=1.0, seed::Int=0, ϵ=1e-14,
		control_amplitude=0.1)
	N = length(Δt)
	speed = zeros(N)
	plin = zeros(N)
	pang = zeros(N)
	ener = zeros(N)
	for i = 1:N
		@show Δt[i]
		speed[i], plin[i], pang[i], ener[i] = mj_batch_astronaut_simulation(; Nsim=Nsim, Δt=Δt[i],
			tsim=tsim, tctrl=tctrl, seed=seed, ϵ=ϵ, control_amplitude=control_amplitude)
	end
    return speed, plin, pang, ener
end


# mj_astronaut_simulation()
# mj_batch_astronaut_simulation(Nsim=10)
Δt = [0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001, 0.00003]
# speed, plin, pang, ener = mj_benchmark_momentum(Δt, Nsim=3)

# Saving results
path = joinpath(module_dir(), "results/astronaut_momentum_mujoco.jld2")
jldsave(path, speed=speed, plin=plin, pang=pang, ener=ener)

plt = plot(layout=(1,2), legend=false, xlims=(1e-2,1e6), ylims=(1e-15,1e6))
plot!(plt[1,1], speed, plin, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
plot!(plt[1,2], speed, pang, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,1], speed, plin, axis=:log, yaxis=:flip, linewidth=3, markersize=6)
scatter!(plt[1,2], speed, pang, axis=:log, yaxis=:flip, linewidth=3, markersize=6)



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
    jd.ctrl .= 1e-5* 0.5 * 10
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
    E[i] = energy(jm, jd)
	Plin[i] = norm(linear_momentum(jm, jd))
	# Pang[i] = angular_momentum(jm, jd)
end
tsim / t
plot(E)
plot(Plin)
plot(Pang)

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel
