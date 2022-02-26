
using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path=dirname(@__DIR__)))
Pkg.instantiate()



using GeometryBasics
using MeshCat
using CoordinateTransformations
using Rotations
using Colors
using LinearAlgebra: I

vis=visualizer()
MeshCat.setobject!(vis[:box1],
    GeometryBasics.Rect3D(Vec(0., 0, 0), Vec(0.1, 0.2, 0.3)))
anim = MeshCat.Animation()
atframe(anim, 0) do
    MeshCat.settransform!(vis[:box1],
        MeshCat.Translation(0., 0, -1) ∘ MeshCat.LinearMap(RotZ(-π/2)))
end
atframe(anim, 30) do
    MeshCat.settransform!(vis[:box1],
        MeshCat.Translation(0., 0, 0) ∘ MeshCat.LinearMap(RotY(π/2)) ∘
            MeshCat.LinearMap(RotZ(π/2)))
end
atframe(anim, 60) do
    MeshCat.settransform!(vis[:box1],
        MeshCat.Translation(0., 0, 1))
end
setanimation!(vis, anim)

open(vis)
render_static(vis)

setobject!(vis[:box1],
    Rect(Vec(0., 0, 0), Vec(0.1, 0.2, 0.3)),


open(joinpath(@__DIR__, "my_scene.html"), "w") do file
    write(file, static_html(vis))
end

MeshPhongMaterial(color=colorant"red"))
