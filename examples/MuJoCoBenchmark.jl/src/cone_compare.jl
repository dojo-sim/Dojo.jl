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

################################################################################
# Demo
################################################################################

Δt = 0.01
jm, jd, mjsim = mj_model("box.xml", Δt=Δt)
tsim = 3.0
N = Int(floor(tsim/Δt))

jd.qpos .= [-1.5,-0.50,0.50,1,0,0,0]
jd.qvel .= [4.0,0.8,0,0,0,0]

traj = zeros(20, 0)
E = zeros(N)
Plin = zeros(N)
Pang = zeros(N)
tcompute = 0.0
ztraj = []
for i = 1:N
    jd.ctrl .= 0.0
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
    E[i] = energy(jm, jd)
	Plin[i] = norm(linear_momentum(jm, jd))
	# Pang[i] = angular_momentum(jm, jd)
	x = deepcopy(jd.qpos[1:3])
	v = deepcopy(jd.qvel[1:3])
	q = deepcopy(jd.qpos[4:7])
	ω = deepcopy(jd.qvel[4:6])
	push!(ztraj, [x; v; q; ω])
end
tsim / t
plot(E)
plot(Plin)
# plot(Pang)

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel

# Saving trajectory
path = joinpath(module_dir(), "results/cone_compare.jld2")
jldsave(path, ztraj=ztraj)
