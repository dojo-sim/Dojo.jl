# Load packages
using MeshCat
using Meshing

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(@__DIR__, "..", "..", "env/tippetop/deps/texture.jl"))

mech = getmechanism(:tippetop, timestep = 0.01, g = -9.00, contact = true, contact_type=:contact);
mech.bodies[3].J = Diagonal([1.9, 2.1, 2.0])
mech.bodies[4].J

initialize!(mech, :tippetop, x = [0,0,1.0], q = UnitQuaternion(RotX(0.01 * π)), ω = [0,0.01,50.])
@elapsed storage = simulate!(mech, 25.0, record = true, verbose = false, opts=SolverOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis=vis)
tippytop_texture!(vis, mech)

setprop!(vis["/Background"], "top_color", RGBA(0.0, 0.0, 0.0))
plot_surface!(vis, x -> x[3] - 0.0)
