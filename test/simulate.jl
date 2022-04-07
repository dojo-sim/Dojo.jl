@testset "step!" begin 
    # get environment and simulate
    env = DojoEnvironments.get_environment("pendulum", 
        timestep=0.1, 
        gravity=0.0);
    reset(env);

    # step (no control)
    z1 = Dojo.get_maximal_state(env.mechanism)
    u1 = zeros(Dojo.input_dimension(env.mechanism))
    z2 = Dojo.step!(env.mechanism, z1, u1)
    @test norm(z2 - z1) < 1.0e-6

    # step (control)
    z1 = Dojo.get_maximal_state(env.mechanism)
    u1 = rand(Dojo.input_dimension(env.mechanism))
    z2 = Dojo.step!(env.mechanism, z1, u1)
    @test norm(z2 - z1) > 1.0e-6
end

@testset "Storage" begin 
    # get environment and simulate
    env = DojoEnvironments.get_environment("pendulum",
        timestep=0.1);
    reset(env);
    DojoEnvironments.initialize_pendulum!(env.mechanism,
        angle=0.25 * π)
    storage = Dojo.simulate!(env.mechanism, 1.0, 
        record=true, 
        verbose=false)
    # check storage length
    @test Dojo.length(storage) == 10
    # convert to vector of vectors
    z = Dojo.get_maximal_state(storage)
    @test length(z) == 10
    # convert back to storage
    s = Dojo.generate_storage(env.mechanism, z)
    @test typeof(s) <: Dojo.Storage
    @test Dojo.length(s) == 10
end
