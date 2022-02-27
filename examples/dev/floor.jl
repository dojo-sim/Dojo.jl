# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis=visualizer()
open(vis)

image = PngImage(joinpath(@__DIR__, "tile_marble_darker.png"))
texture = Texture(image=image, wrap=(1,1), repeat=(14,14))
texture.wrap
texture.repeat
flr_mat = MeshLambertMaterial(map=texture)
flr_mat = MeshPhongMaterial(map=texture)
flr_mat.depthFunc = 3
s = 5.0
flr = HyperRectangle(Vec(0., 0, 0), Vec(s, s, 0.1))
obj = HyperSphere(Point(0.25, 0.25, 0.25), 0.25)
# obj = GeometryBasics.Rect3(Vec(0.0, 0.0), Vec(1.0, 1.0))
obj_mat = MeshPhongMaterial(color=(RGBA(0.75, 0.75, 0.75, 1.0)))
obj_mat.depthFunc =10
setobject!(vis[:floor], flr, flr_mat)
setobject!(vis[:obj], obj, obj_mat)
settransform!(vis[:floor], MeshCat.Translation(-s/2, -s/2, -0.1))




set_floor!(vis, x=5, y=5, tilepermeter=1.0)
setobject!(vis[:obj], obj, obj_mat)
set_light!(vis)



##################################################################################
# Open visualizer
vis=visualizer()
open(vis)



mech = getatlas(timestep=0.01, model_type=:simple)
initializeatlasstance!(mech)
storage = simulate!(mech, 0.20, record=true, verbose=false)
vis, anim = visualize(mech, storage, vis=vis, show_contact=false)

path = joinpath(module_dir(), "env", "atlas", "deps", "mesh", "head_chull.obj")
obj = MeshFileObject(path)
# obj = HyperRectangle(Vec(0.,0.,0.), Vec(1.,1.,1.))
mat = MeshPhongMaterial(color=RGBA(1,1,1,1))
setobject!(vis[:head0], obj)
settransform!(vis[:head0], MeshCat.Translation(0.0, 0.5, 0.0))

image = PngImage(joinpath(module_dir(), "assets", "tile.png"))
texture = Texture(image=image, wrap=(1,1), repeat=(1,1))
mat = MeshPhongMaterial(map=texture)
obj = HyperRectangle(Vec(0., 0, 0), Vec(1, 1, 1))
setobject!(vis[:floor], obj, mat)
setobject!(vis[:floor], obj, mat)
settransform!(vis[:floor], MeshCat.Translation(-x/2, -y/2, -z))

vis
