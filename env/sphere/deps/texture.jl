function sphere_texture!(vis::Visualizer, mech::Mechanism;
        imagepath = joinpath(@__DIR__, "dojo_texture.png"), name=:robot)
    image = PngImage(imagepath)
    texture = Texture(image=image);
    material = MeshLambertMaterial(map=texture);
    setobject!(vis[name]["bodies"]["body:1"], convert_shape(mech.bodies.values[1].shape), material)
    return nothing
end
