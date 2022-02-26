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

timestep=0.01
jm, jd, mjsim = mj_model("dzhanibekov.xml", timestep=timestep)
tsim = 1.0
N = Int(floor(tsim/timestep))

# jd.qpos
jd.qvel .= [0,0,0,20π,0,0.01]

traj = zeros(20, 0)
E = zeros(N)
Plin = zeros(N)
Pang = zeros(N)
tcompute = 0.0
for i = 1:N
    jd.ctrl .= 0.0
    tcompute += @elapsed mj_step(jm, jd);
    traj = hcat(traj, deepcopy(getstate(mjsim)))
    E[i] = energy(jm, jd)
	Plin[i] = norm(linear_momentum(jm, jd))
	# Pang[i] = angular_momentum(jm, jd)
end
tsim / tcompute
plot(E)
plot(Plin)
# plot(Pang)

LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)

jd.qpos
jd.qvel


timestep=0.01
jm, jd, mjsim = mj_model("dzhanibekov.xml", timestep=timestep)
jd.qvel .= [0,0,0,20π,0,0]
# jd.qvel .= [1,0,0,0.0,0,0]

mj_step(jm, jd)
MuJoCo.mj_subtreeVel(jm, jd)
jd.subtree_linvel
jd.subtree_angmom

jm.body_mass
mjsim.d.xipos

function linear_momentum(jm, jd)
	MuJoCo.mj_subtreeVel(jm, jd)
    p = Vector(sum(jd.subtree_linvel, dims=2)[:,1])
    return p
end






#
# ################################################################################
# ################################################################################
#
# timestep=0.01
# jm, jd, mjsim = mj_model("half_cheetah.xml", timestep=timestep)
# jd.qpos[2] += 2.0
# jd.qvel .= [0,0,0, 0,0,0, 0,0,5]
# # mj_step(jm, jd)
# MuJoCo.mj_subtreeVel(jm, jd)
# jd.subtree_linvel
# jd.subtree_angmom[:,1:4]
# jd.subtree_angmom[:,5:8]
# traj = zeros(28, 0)
#
#
#
# #
# # jm, jd, mjsim = mj_model("astronaut.xml", timestep=timestep)
# # tsim = 0.10
# # N = Int(floor(tsim/timestep))
#
#
# # jd.qpos
# # jd.qvel .= [1; zeros(26)]
#
# # traj = zeros(83, 0)
# E = zeros(N)
# Plin = zeros(N)
# Pang = zeros(N)
# tcompute = 0.0
# for i = 1:N
#     jd.ctrl .= 0.0
#     tcompute += @elapsed mj_step(jm, jd);
#     traj = hcat(traj, deepcopy(getstate(mjsim)))
#     E[i] = energy(jm, jd)
# 	Plin[i] = norm(linear_momentum(jm, jd))
# 	# Pang[i] = angular_momentum(jm, jd)
# end
# tsim / tcompute
# plot(E)
# plot(Plin)
# # plot(Pang)
#
# LyceumMuJoCoViz.visualize(mjsim, trajectories=traj)
#
#
#
# MuJoCo.mj_subtreeVel(jm, jd)
# jd.subtree_linvel
# jd.subtree_angmom
#
#
