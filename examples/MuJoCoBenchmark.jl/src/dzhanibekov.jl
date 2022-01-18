# Utils
function module_dir()
    return joinpath(@__DIR__, "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using LinearAlgebra
using MuJoCo
using LyceumBase
using LyceumMuJoCo
using LyceumMuJoCoViz
using Plots

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
m = jlModel(joinpath(module_dir(), "model/tennis_racket.xml"))
dt = 0.025
N = 15000
H = N*dt/60
MuJoCo.@set!! m.opt.timestep = dt
d = jlData(m)
mj = MJSim(m, d)


d.qvel .= [0,0,0,20π,0,0.00]
# d.qvel .= [0.5,0,0,0,0,0]
traj = zeros(20,0)
E = []
for i=1:N
    mj_step(m, d);
    traj = hcat(traj, deepcopy(getstate(mj)))
    push!(E, energy(m, d.qpos, d.qvel))
end
plot(E)

LyceumMuJoCoViz.visualize(mj, trajectories=traj)
