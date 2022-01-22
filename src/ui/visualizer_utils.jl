
function set_floor!(vis::Visualizer; x=20.0, y=20.0, z=0.1, color=nothing,
        tilepermeter=1.0, imagename="tile.png")
    image = PngImage(joinpath(module_dir(), "assets", imagename))
    repeat = Int.(ceil.(tilepermeter * [x,y]))
    texture = Texture(image=image, wrap=(1,1), repeat=(repeat[1],repeat[2]))
    mat = MeshPhongMaterial(map=texture)
    (color != nothing) && (mat = MeshPhongMaterial(color=color))
    obj = HyperRectangle(Vec(0., 0, 0), Vec(x, y, z))
    setobject!(vis[:floor], obj, mat)
    settransform!(vis[:floor], MeshCat.Translation(-x/2, -y/2, -z))
    return nothing
end

function set_surface!(vis::Visualizer, f::Any; xlims = [-20.0, 20.0],
    ylims = [-20.0, 20.0], color=RGBA(1.0, 1.0, 1.0, 0.0), n::Int=200)
    mesh = GeometryBasics.Mesh(f,
        HyperRectangle(Vec(xlims[1], ylims[1], -2.0), Vec(xlims[2]-xlims[1], ylims[2]-ylims[1], 4.0)),
        Meshing.MarchingCubes(), samples=(n, n, Int(floor(n/8))))
    setobject!(vis["surface"], mesh, MeshPhongMaterial(color=color))
    return nothing
end

function set_light!(vis::Visualizer; ambient=0.35, fill=0.25, pointX=0.85,
        pointXshadow::Bool=true, direction::String="Positive")
    setprop!(vis["/Lights/AmbientLight/<object>"], "intensity", ambient)
    setprop!(vis["/Lights/FillLight/<object>"], "intensity", 0.25)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "intensity", 0.85)
    setprop!(vis["/Lights/PointLight$(direction)X/<object>"], "castShadow", true)
    return nothing
end
