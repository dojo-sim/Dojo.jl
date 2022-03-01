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
# Demo Elliptic
################################################################################

timestep = 0.01
jm, jd, mjsim = mj_model("box_elliptic.xml", timestep=timestep)
tsim = 3.0
N = Int(floor(tsim/timestep))

jd.qpos .= [-1.5,-0.50,0.50,1,0,0,0]
jd.qvel .= [4.0,0.8,0,0,0,0]

traj = zeros(20, 0)
tcompute = 0.0
ztraj = []
for i = 1:N
    jd.ctrl .= 0.0
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
	x = deepcopy(jd.qpos[1:3])
	v = deepcopy(jd.qvel[1:3])
	q = deepcopy(jd.qpos[4:7])
	ω = deepcopy(jd.qvel[4:6])
	push!(ztraj, [x; v; q; ω])
end
tsim / tcompute

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel

# Saving trajectory
path = joinpath(module_dir(), "results/cone_compare_elliptic.jld2")
jldsave(path, ztraj=ztraj)

################################################################################
# Demo Pyramidal
################################################################################

timestep = 0.01
jm, jd, mjsim = mj_model("box_pyramidal.xml", timestep=timestep)
tsim = 3.0
N = Int(floor(tsim/timestep))

jd.qpos .= [-1.5,-0.50,0.50,1,0,0,0]
jd.qvel .= [4.0,0.8,0,0,0,0]

traj = zeros(20, 0)
tcompute = 0.0
ztraj = []
for i = 1:N
    jd.ctrl .= 0.0
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
	x = deepcopy(jd.qpos[1:3])
	v = deepcopy(jd.qvel[1:3])
	q = deepcopy(jd.qpos[4:7])
	ω = deepcopy(jd.qvel[4:6])
	push!(ztraj, [x; v; q; ω])
end
tsim / tcompute

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel

# Saving trajectory
path = joinpath(module_dir(), "results/cone_compare_pyramidal.jld2")
jldsave(path, ztraj=ztraj)
ztraj
