@testset "Pendulum" begin
    # get pendulum environment and simulate
    timestep=0.1
    env = DojoEnvironments.get_environment("pendulum", 
        timestep=timestep, 
        gravity=-10.0);
    reset(env);
    DojoEnvironments.initialize_pendulum!(env.mechanism, 
        angle=0.25 * π)
    u1 = rand(Dojo.input_dimension(env.mechanism))
    z1 = Dojo.get_maximal_state(env.mechanism)

    storage = Dojo.simulate!(env.mechanism, 1.0, 
        record=true)
    @test norm(z1 - Dojo.get_maximal_state(storage, 1)) < 1.0e-6
    zTs = Dojo.get_maximal_state(env.mechanism)
    @test norm(zTs - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6
    zT = Dojo.get_next_state(env.mechanism)
    DojoEnvironments.initialize_pendulum!(env.mechanism, 
        angle=0.25 * π)
    storage = Dojo.simulate!(env.mechanism, 1.0 + timestep, 
        record=true)
    @test norm(zT - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6

    # maximal gradients
    Dojo.step!(env.mechanism, z1, u1)
    Fx1, Fu1 = Dojo.get_maximal_gradients(env.mechanism)
    Fx2, Fu2 = Dojo.get_maximal_gradients!(env.mechanism, z1, u1)
    @test norm(Fx1 - Fx2, Inf) < 1.0e-6 
    @test norm(Fu1 - Fu2, Inf) < 1.0e-6 

    timestep=0.1
    env = DojoEnvironments.get_environment("halfcheetah", 
        timestep=timestep, 
        gravity=-10.0);
    reset(env);

    # get body 
    @test Dojo.get_body(env.mechanism, :origin).name == :origin

    # get contact
    contact_name = env.mechanism.contacts[1].name
    @test Dojo.get_contact(env.mechanism, contact_name).name == contact_name

    # velocity indices 
    @test Dojo.velocity_index(env.mechanism) == [4, 5, 6, 8, 10, 12, 14, 16, 18]

    # set state
    z = Dojo.get_maximal_state(env.mechanism)
    env2 = DojoEnvironments.get_environment(:halfcheetah, 
        timestep=timestep, 
        gravity=-10.0);
    reset(env);
    Dojo.set_maximal_state!(env2.mechanism, zeros(Dojo.maximal_dimension(env2.mechanism)))
    z2 = Dojo.get_maximal_state(env2.mechanism)
    @test norm(z2) < 1.0e-8
    Dojo.set_maximal_state!(env2.mechanism, z)
    @test norm(z - Dojo.get_maximal_state(env2.mechanism)) < 1.0e-8

    # change atlas floating base 
    mech = get_mechanism(:atlas)
    j = get_joint(mech, :floating_base)
    base_origin = get_body(mech, j.child_id).name 

    mech = set_floating_base(mech, :head) 
    j = get_joint(mech, :floating_base)
    base_new = get_body(mech, j.child_id).name 

    @test base_origin != base_new 
    @test base_new == :head
end

#TODO: get and set methods
