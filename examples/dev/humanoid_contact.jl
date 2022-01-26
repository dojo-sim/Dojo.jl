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

gravity = -1.0
dt = 0.01
cf = 0.8
damper = 4.5
spring = 0.0
mech = gethumanoid(g=gravity, timestep=dt, cf=cf, damper=damper, spring=spring)
initialize!(mech, :humanoid, rot = [0,0,0.3], tran = [0,0,1.4])
storage = simulate!(mech, 1.0, record = true, verbose = true,
	opts = SolverOptions(rtol = 1e-10, btol = 1e-6, verbose = false))
visualize(mech, storage, vis=vis, show_contact = true)
mech.bodies[14]

# v = 14
# c = 27
# mech.bodies[v]
# mech.contacts[c]
# matrix_entries = mech.system.matrix_entries
# diagonal_inverses = mech.system.diagonal_inverses
# matrix_entries[v,v]
# matrix_entries[v,c]
# matrix_entries[c,c]
# matrix_entries[c,v]
# diagonal_inverses[c]
# inv(matrix_entries[c,c].value)
# rank(matrix_entries[c,c].value)



vis = Visualizer()
open(vis)

object = Cylinder(0.2, 0.3)
setobject!(vis, object)

setobject!(subvisshape, visshape, MeshPhongMaterial(color=(transparent ? RGBA(0.75, 0.75, 0.75, 0.5) : shape.color)))
