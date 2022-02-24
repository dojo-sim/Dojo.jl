@testset "Bodies: Shapes" begin
    # shapes
    box = Box(1.0, 1.0, 1.0)
    cylinder = Cylinder(1.0, 1.0)
    capsule = Capsule(1.0, 1.0)
    shapes = Shapes([box, cylinder, capsule])
    sphere = Sphere(1.0)
    pyramid = Pyramid(1.0, 1.0)
    mesh = Mesh(joinpath(@__DIR__, "../env/atlas/deps/mesh/head.obj"))

    # convert 
    box_geom = convert_shape(box)
    cylinder_geom = convert_shape(cylinder)
    capsule_geom = convert_shape(capsule)
    shapes_geom = convert_shape(shapes)
    sphere_geom = convert_shape(sphere)
    pyramid_geom = convert_shape(pyramid)
    mesh_geom = convert_shape(mesh)
    @test true
end


