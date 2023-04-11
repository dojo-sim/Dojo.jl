@testset "Utilities" begin
    # create visualizer
    vis = Dojo.Visualizer();

    # default background 
    Dojo.set_background!(vis)

    # add floor 
    Dojo.set_floor!(vis)

    # add quadratic surface 
    Dojo.set_surface!(vis, x -> x' * x, 
        color=Dojo.RGBA(1.0, 0.0, 0.0, 1.0))

    # lighting 
    Dojo.set_light!(vis)

    # camera 
    Dojo.set_camera!(vis, 
        zoom=0.25)

    # test that methods don't fail
    @test true
end

@testset "MeshCat mechanism" begin
    # create visualizer
    vis = Dojo.Visualizer();
    # get mechanism and simulate
    mechanism = DojoEnvironments.get_mechanism(:halfcheetah;
        timestep=0.1)
    storage = Dojo.simulate!(mechanism, 0.25;
        record=true, verbose=true);
    # visualize simulation
    Dojo.visualize(mechanism, storage; vis=vis)

    # test that methods don't fail
    @test true
end

@testset "URDF mesh" begin
    # create visualizer
    vis = Dojo.Visualizer();
    # get mechanism and simulate
    mechanism = DojoEnvironments.get_mechanism(:quadruped; timestep=0.1)
    storage = Dojo.simulate!(mechanism, 0.25;
        record=true, verbose=false)
    # visualize simulation
    Dojo.visualize(mechanism, storage; vis=vis)
    # visualize w/ contact
    Dojo.visualize(mechanism, storage; 
        vis=vis, show_contact=true)

    # test that methods don't fail
    @test true
end
