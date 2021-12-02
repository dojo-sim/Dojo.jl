# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat
# using IterativeLQR

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))



mech = getmechanism(:quadruped, Δt = Δt, g = gravity, cf = cf, damper = 0*10.0, spring = 0*300.0)
initialize!(mech, :quadruped)
@elapsed storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
