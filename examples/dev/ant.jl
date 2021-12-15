# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


mech = getmechanism(:ant, Δt = 0.05, g = -9.81, contact = true,
    contact_body = true, spring = 0.0, damper = 1.0);
initialize!(mech, :ant, rot = [0,0,0.], ankle = 0.25)
@elapsed storage = simulate!(mech, 2.0, record = true, verbose = false,
    opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)

env = make("ant", vis = vis)

env.aspace
seed(env, s = 11)
obs = reset(env)[2]
render(env)

1000*sample(env.aspace)
collect(env.mechanism.eqconstraints)[1]
for i = 1:25
    render(env)
    sleep(0.05)
    # action = 120*env.mechanism.Δt*ones(6)#1000*sample(env.aspace) # your agent here (this takes random actions)
    action = sample(env.aspace)#1000*sample(env.aspace) # your agent here (this takes random actions)
    obs, r, done, info = step(env, action)
    @show r

    if done
        observation = reset(env)
    end
end
close(env)

env.mechanism.eqconstraints
controldim(env.mechanism)
sample(env.aspace)
# sample(env.aspace)
