# Load packages
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(@__DIR__, "..", "..", "env/tippetop/deps/texture.jl"))

mech = getmechanism(:tippetop, Δt = 0.01, g = -9.00, contact = true, contact_mode=:soc);
mech.bodies[3].J = Diagonal([1.9, 2.1, 2.0])
mech.bodies[4].J

initialize!(mech, :tippetop, x = [0,0,1.0], q = UnitQuaternion(RotX(0.01 * π)), ω = [0,0.01,50.])
@elapsed storage = simulate!(mech, 25.0, record = true, verbose = false, opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis=vis)
tippytop_texture!(vis, mech)

setprop!(vis["/Background"], "top_color", RGBA(0.0, 0.0, 0.0))

function plot_surface!(vis::Visualizer, f::Any; xlims = [-10.0, 10.0],
    ylims = [-10.0, 10.0], col=[1.0, 215/255, 0.0], α=1.0, n::Int=200)
    mesh = GeometryBasics.Mesh(f,
        HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2]-xlims[1], ylims[2]-ylims[1], 4.0)),
        Meshing.MarchingCubes(), samples=(n, n, Int(floor(n/8))))
    setobject!(vis["surface"], mesh,
            MeshPhongMaterial(color=RGBA{Float32}(col..., α)))
    return nothing
end

using Meshing
plot_surface!(vis, x -> x[3] - 0.0)

mech = getmechanism(:rectangle, Δt = 0.01, g = -0.00);
body1 = collect(mech.bodies)[1]
body1.J = Matrix(Diagonal([1,2,3.]))
initialize!(mech, :rectangle, x = [0,0,1.0], q = UnitQuaternion(RotX(0.00*π)), ω = [0,5.0,0.01])
@elapsed storage = simulate!(mech, 8.0, record = true, verbose = false, opts=InteriorPointOptions(verbose=false, btol = 1e-6))
visualize(mech, storage, vis = vis)
