function tippytop_texture!(vis::Visualizer, mech::Mechanism; 
    image1=joinpath(@__DIR__, "escher.png"), image2=joinpath(@__DIR__, "escher.png"))
    image = PngImage(image1)
    texture = Texture(image=image);   
    material = MeshLambertMaterial(map=texture);       
    setobject!(vis["bodies"]["body:1"], convert_shape(mech.bodies[3].shape), material)
    
    image = PngImage(image2)
    texture = Texture(image=image);   
    material = MeshLambertMaterial(map=texture); 
    
    setobject!(vis["bodies"]["body:2"], convert_shape(mech.bodies[4].shape), material)
end
