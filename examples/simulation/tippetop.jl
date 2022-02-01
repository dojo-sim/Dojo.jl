using Dojo
using LinearAlgebra 

# ## Open visualizer
vis = Visualizer()
open(vis)

# ## Include model
include(joinpath(@__DIR__, "..", "..", "env/tippetop/deps/texture.jl"))
include(joinpath(@__DIR__, "..", "..", "env/tippetop/methods/initialize.jl"))

mech = get_mechanism(:tippetop, timestep = 0.01, gravity=-9.81, contact = true, contact_type=:contact)
initialize!(mech, :tippetop, x = [0,0,1.0], q = UnitQuaternion(RotX(0.01 * π)), ω = [0.0, 0.01, 50.0])
@elapsed storage = simulate!(mech, 25.0, record = true, verbose = false, opts=SolverOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis=vis)

# tippytop_texture!(vis, mech)
# setprop!(vis["/Background"], "top_color", RGBA(0.0, 0.0, 0.0))
# plot_surface!(vis, x -> x[3] - 0.0)
with