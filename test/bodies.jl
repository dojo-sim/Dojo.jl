@testset "Bodies: Shapes" begin
    # shapes
    box = Dojo.Box(1.0, 1.0, 1.0)
    cylinder = Dojo.Cylinder(1.0, 1.0)
    capsule = Dojo.Capsule(1.0, 1.0)
    shapes = Dojo.Shapes([box, cylinder, capsule])
    sphere = Dojo.Sphere(1.0)
    pyramid = Dojo.Pyramid(1.0, 1.0)
    mesh = Dojo.Mesh(joinpath(@__DIR__, "../environments/atlas/deps/mesh/head.obj"))

    # convert 
    box_geom = Dojo.convert_shape(box)
    cylinder_geom = Dojo.convert_shape(cylinder)
    capsule_geom = Dojo.convert_shape(capsule)
    shapes_geom = Dojo.convert_shape(shapes)
    sphere_geom = Dojo.convert_shape(sphere)
    pyramid_geom = Dojo.convert_shape(pyramid)
    mesh_geom = Dojo.convert_shape(mesh)
    @test true
end


