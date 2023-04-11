@testset "step!" begin 
    # get mechanism and simulate
    mechanism = DojoEnvironments.get_mechanism(:pendulum;
        timestep=0.1, 
        gravity=0.0);

    # step (no control)
    z1 = Dojo.get_maximal_state(mechanism)
    u1 = zeros(Dojo.input_dimension(mechanism))
    z2 = Dojo.step!(mechanism, z1, u1)
    @test norm(z2 - z1) < 1.0e-6

    # step (control)
    z1 = Dojo.get_maximal_state(mechanism)
    u1 = rand(Dojo.input_dimension(mechanism))
    z2 = Dojo.step!(mechanism, z1, u1)
    @test norm(z2 - z1) > 1.0e-6
end

@testset "Storage" begin 
    # get mechanism and simulate
    mechanism = DojoEnvironments.get_mechanism(:pendulum;
        timestep=0.1);
    DojoEnvironments.initialize_pendulum!(mechanism,
        angle=0.25 * π)
    storage = Dojo.simulate!(mechanism, 1.0, 
        record=true, 
        verbose=false)
    # check storage length
    @test Dojo.length(storage) == 10
    # convert to vector of vectors
    z = Dojo.get_maximal_state(storage)
    @test length(z) == 10
    # convert back to storage
    s = Dojo.generate_storage(mechanism, z)
    @test typeof(s) <: Dojo.Storage
    @test Dojo.length(s) == 10
end
