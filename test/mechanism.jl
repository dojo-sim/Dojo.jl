@testset "Pendulum" begin
    # get pendulum mechanism and simulate
    timestep=0.1
    mechanism = DojoEnvironments.get_mechanism(:pendulum;
        timestep, gravity=-10.0);
    DojoEnvironments.initialize!(mechanism, :pendulum; angle=0.25 * π)
    u1 = rand(Dojo.input_dimension(mechanism))
    z1 = Dojo.get_maximal_state(mechanism)

    storage = Dojo.simulate!(mechanism, 1.0; record=true)
    @test norm(z1 - Dojo.get_maximal_state(storage, 1)) < 1.0e-6
    zTs = Dojo.get_maximal_state(mechanism)
    @test norm(zTs - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6
    zT = Dojo.get_next_state(mechanism)
    DojoEnvironments.initialize!(mechanism, :pendulum; angle=0.25 * π)
    storage = Dojo.simulate!(mechanism, 1.0 + timestep; record=true)
    @test norm(zT - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6

    # maximal gradients
    Dojo.step!(mechanism, z1, u1)
    Fx1, Fu1 = Dojo.get_maximal_gradients(mechanism)
    Fx2, Fu2 = Dojo.get_maximal_gradients!(mechanism, z1, u1)
    @test norm(Fx1 - Fx2, Inf) < 1.0e-6 
    @test norm(Fu1 - Fu2, Inf) < 1.0e-6 

    timestep=0.1
    mechanism = DojoEnvironments.get_mechanism(:halfcheetah;
        timestep, gravity=-10.0);

    # get body 
    @test Dojo.get_body(mechanism, :origin).name == :base

    # get contact
    contact_name = mechanism.contacts[1].name
    @test Dojo.get_contact(mechanism, contact_name).name == contact_name

    # velocity indices 
    @test Dojo.velocity_index(mechanism) == [4, 5, 6, 8, 10, 12, 14, 16, 18]

    # set state
    z = Dojo.get_maximal_state(mechanism)
    mechanism = DojoEnvironments.get_mechanism(:halfcheetah;
        timestep, gravity=-10.0);
    Dojo.set_maximal_state!(mechanism, zeros(Dojo.maximal_dimension(mechanism)))
    z2 = Dojo.get_maximal_state(mechanism)
    @test norm(z2) < 1.0e-8
    Dojo.set_maximal_state!(mechanism, z)
    @test norm(z - Dojo.get_maximal_state(mechanism)) < 1.0e-8

    # change atlas floating base 
    mechanism = get_mechanism(:atlas; parse_dampers=false)
    j = get_joint(mechanism, :floating_base)
    base_origin = get_body(mechanism, j.child_id).name 

    mechanism = set_floating_base(mechanism, :head) 
    j = get_joint(mechanism, :floating_base)
    base_new = get_body(mechanism, j.child_id).name 

    @test base_origin != base_new 
    @test base_new == :head
end

@testset "Remove fixed joints" begin
    origin = Origin()
    box1 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box2 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box3 = Box(0.1,0.1,0.5,0.5;color=RGBA(1,1,0,0.5))
    box4 = Box(0.1,0.1,0.5,0.5;color=RGBA(1,0,1,0.5))
    joint1 = JointConstraint(Revolute(origin, box1, [1;0;0], child_vertex = [0;0;0.25]))
    joint2 = JointConstraint(Fixed(box1, box2; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    joint3 = JointConstraint(Revolute(origin, box3, [1;0;0], child_vertex = [0;0;0.25]))
    joint4 = JointConstraint(Fixed(box3, box4; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    mechanism = Mechanism(origin, [box1;box2;box3;box4], [joint1;joint2;joint3;joint4]) 

    origin, bodies, joints = Dojo.reduce_fixed_joints(mechanism.origin,mechanism.bodies,mechanism.joints)
    box1 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box2 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    joint1 = JointConstraint(Revolute(origin, box1, [1;0;0]; child_vertex = [0;0;0.25]))
    joint2 = JointConstraint(Fixed(box1, box2; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    mechanism = Mechanism(origin, [box1,box2,bodies[2]], [joint1,joint2,joints[2]]; gravity = -9.81, timestep = 0.01)

    Dojo.zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism,joint1,[pi/2])
    set_minimal_coordinates!(mechanism,joints[2],[pi/2])
    
    storage = simulate!(mechanism, 5.0; record = true)

    @test abs(storage.q[3][end].v1 - storage.q[1][end].v1) < 1.0e-5
end

@testset "Force and impulse input" begin
    timestep = 0.01

    # Force input, default scaling
    origin0 = Origin()
    box0 = Box(1.0,1,1,1)
    joint0 = JointConstraint(Floating(origin0,box0))
    mech0 = Mechanism(origin0, [box0], [joint0]; gravity = -9.81, timestep)
    controller0!(mechanism, k) = set_input!(mechanism, [0;0;9.81;0;0;0])
    storage0 = simulate!(mech0, 10.0, controller0!; record = true)

    # Force input, default scaling explicitly set
    origin1 = Origin()
    box1 = Box(1.0,1,1,1)
    joint1 = JointConstraint(Floating(origin1,box1))
    mech1 = Mechanism(origin1, [box1], [joint1]; gravity = -9.81, timestep, input_scaling = timestep)
    controller1!(mechanism, k) = set_input!(mechanism, [0;0;9.81;0;0;0])
    storage1 = simulate!(mech1, 10.0, controller1!; record = true)

    # Impulse input, input scaled by timestep
    origin2 = Origin()
    box2 = Box(1.0,1,1,1)
    joint2 = JointConstraint(Floating(origin2,box2))
    mech2 = Mechanism(origin2, [box2], [joint2]; gravity = -9.81, timestep, input_scaling = 1)
    controller2!(mechanism, k) = set_input!(mechanism, [0;0;9.81*timestep;0;0;0])
    storage2 = simulate!(mech2, 10.0, controller2!; record = true)

    @test norm(storage0.x[1][end] - storage1.x[1][end]) < 1.0e-5
    @test norm(storage1.x[1][end] - storage2.x[1][end]) < 1.0e-5
end


#TODO: get and set methods
