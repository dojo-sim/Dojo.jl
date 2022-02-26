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

timestep = 0.01
jm, jd, mjsim = mj_model("humanoid.xml", timestep=0.01)
tsim = 1.0
N = Int(floor(tsim/timestep))

jd.qpos[3] += 1.0
Random.seed!(100)
jd.qvel .= 1.0 * (rand(27) .- 0.5)

traj = zeros(83, 0)
E = zeros(N)
Plin = zeros(N)
Pang = zeros(N)
tcompute = 0.0
for i = 1:N
    jd.ctrl .= 1e-5* 0.5
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
