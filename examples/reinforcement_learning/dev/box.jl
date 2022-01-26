################################################################################
# Script
################################################################################

env = make("dice", timestep = 0.05, g = -9.81);
reset(env, x = [0,0,0.2])
action = nothing
step(env, action)

env = make("dice");
for i = 1:20
    observation = reset(env)
    for t = 1:200
        # render(env)
        println(scn.(observation))
        action = sample(env.aspace)
        observation, reward, done, info = step(env, action)
        if done
            println("Episode finished after $(t+1) timesteps")
            break
        end
    end
end
close(env)

################################################################################
# Development
################################################################################
# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

mech = getmechanism(:dice, timestep = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, contact_type = :contact)
# mech = getmechanism(:dice, timestep = 0.01, g = -9.81, cf = 0.2, contact = true, mode=:box, contact_type = :linear_contact)
# mech = getmechanism(:dice, timestep = 0.01, g = -9.81, contact = true, mode=:box, contact_type = :impact)
Random.seed!(100)
ω = 0.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = 1.0 * [4, 1, 0.0]
initialize!(mech, :dice, x = x, v = v, ω = ω)
storage = simulate!(mech, 1.3, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
