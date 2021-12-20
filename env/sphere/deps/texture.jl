function sphere_texture!(vis::Visualizer, mech::Mechanism;
        image1 = joinpath(@__DIR__, "dojo_texture.png"))
    image = PngImage(image1)
    texture = Texture(image=image);
    material = MeshLambertMaterial(map=texture);
    setobject!(vis["bodies"]["body:1"], convertshape(mech.bodies.values[1].shape), material)
    return nothing
end
