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

function energy(m, qpos, qvel)
    M = m.body_mass[2]
    J = Diagonal(m.body_inertia[:,2])
    x = qpos[1:3]
    q = qpos[4:7]
    v = qvel[1:3]
    ω = qvel[4:6]
    E = 0.5*M*v'*v
    E += 0.5*ω'*J*ω
    return E
end

module_dir()
m = jlModel(joinpath(module_dir(), "model/box.xml"))
dt = 0.01
N = 300
H = N*dt
MuJoCo.@set!! m.opt.timestep = dt
d = jlData(m)
mj = MJSim(m, d)

d.qpos
d.qpos .= [-1.5,-0.50,0.50,1,0,0,0]
d.qvel .= [4.0,0.8,0,0,0,0]
traj = zeros(20,0)
ztraj = []
E = []
for i=1:N
    mj_step(m, d);
    traj = hcat(traj, deepcopy(getstate(mj)))
    x = deepcopy(d.qpos[1:3])
    v = deepcopy(d.qvel[1:3])
    q = deepcopy(d.qpos[4:7])
    ω = deepcopy(d.qvel[4:6])
    push!(ztraj, [x; v; q; ω])
    push!(E, energy(m, d.qpos, d.qvel))
end
plot(E)

LyceumMuJoCoViz.visualize(mj, trajectories=traj)

d.qpos
d.qvel

# Saving trajectory
path = joinpath(module_dir(), "results/cone_compare.jld2")
jldsave(path, ztraj=ztraj)
