# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

################################################################################
# Dzhanibekov effect
################################################################################
Δt0 = 0.01
g0 = 0.0
mech = getmechanism(:dzhanibekov, Δt = Δt0, g = g0)

bodies = collect(mech.bodies)
bodies[1].m
bodies[2].m
bodies[1].J
bodies[2].J

initialize!(mech, :dzhanibekov, x = [0,0,0.], v = [0,0,0.], q = UnitQuaternion(1.,0,0,0), ω = [5,0.5,0.])

storage = simulate!(mech, 10.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)


################################################################################
# Full Friction Cone vs. Linearized Friction Cone
################################################################################
Δt0 = 0.01
g0 = -9.81
cf0 = 0.30
mech = getmechanism(:dice, Δt = Δt0, g = g0, cf = cf0)

bodies = collect(mech.bodies)
bodies[1].m
bodies[1].J

initialize!(mech, :dice, x = [0,0,1.], v = [4,2,1.], q = UnitQuaternion(1.,0,0,0), ω = [0,0,0.])

storage = simulate!(mech, 10.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)
