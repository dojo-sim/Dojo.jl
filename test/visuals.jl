@testset "Visualizer: Utilities" begin
    # create visualizer
    vis = Dojo.Visualizer();

    # default background 
    Dojo.set_background!(vis)

    # add floor 
    Dojo.set_floor!(vis)

    # add quadratic surface 
    Dojo.set_surface!(vis, x -> x' * x, color=Dojo.RGBA(1.0, 0.0, 0.0, 1.0))

    # lighting 
    Dojo.set_light!(vis)

    # camera 
    Dojo.set_camera!(vis, zoom=0.25)

    # test that methods don't fail
    @test true
end

@testset "Visualizer: MeshCat mechanism" begin
    # create visualizer
    vis = Dojo.Visualizer();
    # get environment and simulate
    env = Dojo.make("halfcheetah", dt=0.1)
    Dojo.reset(env)
    storage = Dojo.simulate!(env.mechanism, 0.25, record=true, verbose=true);
    # visualize simulation
    Dojo.visualize(env.mechanism, storage, vis=vis)
end

@testset "Visualizer: URDF mesh" begin
    # create visualizer
    vis = Dojo.Visualizer();
    # get environment and simulate
    env = Dojo.make("quadruped", dt=0.1)
    Dojo.reset(env)
    storage = Dojo.simulate!(env.mechanism, 0.25, record=true, verbose=false)
    # visualize simulation
    Dojo.visualize(env.mechanism, storage, vis=vis)
    # visualize w/ contact
    Dojo.visualize(env.mechanism, storage, vis=vis, show_contact=true)
    # test that methods don't fail
    @test true
end
