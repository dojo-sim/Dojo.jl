using Pkg; Pkg.activate(@__DIR__)
# add Lyceum registry: add registry https://github.com/Lyceum/LyceumRegistry.git
using MuJoCo
mj_activate("/home/taylor/.mujoco/bin/mjkey.txt") # set location to MuJoCo key path

using LyceumMuJoCo, LyceumMuJoCoViz 

# ## load MuJoCo 
path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5_10Hz.xml")
# path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5_50Hz.xml")
# path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5_100Hz.xml")
# path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5_250Hz.xml")
# path = joinpath(@__DIR__, "../../../env/atlas/deps/atlas_v5_1000Hz.xml")

include("mujoco_model.jl")
atlas = MuJoCoModel(path)
sim = LyceumMuJoCo.MJSim(atlas.m, atlas.d)

T = 10 * 1
states = Array(undef, statespace(sim), T)
# sim.d.qpos = copy(sim.m.qpos0)
sim.d.qpos[3] += 0.5 # = zero(sim.d.qpos)
@show sim.d.qpos
sim.d.qvel .= 0.0
max_viol = Float64[]
for t = 1:T
    LyceumMuJoCo.step!(sim)
    states[:, t] .= getstate(sim)
    push!(max_viol, !isempty(sim.d.contact) ? minimum([c.dist for c in sim.d.contact]) : 0.0) 
end

visualize(sim, trajectories=[states])

@show minimum(max_viol)

