@testset "Shape convertion" begin
    # shapes
    box = Dojo.Box(1.0, 1.0, 1.0)
    cylinder = Dojo.Cylinder(1.0, 1.0)
    capsule = Dojo.Capsule(1.0, 1.0)
    shapes = Dojo.CombinedShapes([box, cylinder, capsule])
    sphere = Dojo.Sphere(1.0)
    pyramid = Dojo.Pyramid(1.0, 1.0)
    frame = Dojo.FrameShape()
    mesh = Dojo.Mesh(joinpath(dirname(pathof(DojoEnvironments)), "../mechanisms/atlas/dependencies/mesh/head.obj"))

    # convert
    box_geom = Dojo.convert_shape(box)
    cylinder_geom = Dojo.convert_shape(cylinder)
    capsule_geom = Dojo.convert_shape(capsule)
    shapes_geom = Dojo.convert_shape(shapes)
    sphere_geom = Dojo.convert_shape(sphere)
    pyramid_geom = Dojo.convert_shape(pyramid)
    frame_geom = Dojo.convert_shape(frame)
    mesh_geom = Dojo.convert_shape(mesh)
    @test true
end
