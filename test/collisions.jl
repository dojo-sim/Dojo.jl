@testset "Collision: Sphere-sphere (nonlinear)" begin
    # Parameters

    function get_two_body(; 
        radius_body1=0.5, 
        radius_body2=0.5,
        mass_body1=1.0, 
        mass_body2=1.0,
        friction_type=:nonlinear,
        friction_coefficient=0.5,
        joint_world_body1=:Fixed,
        gravity=-9.81,
        timestep=0.1)

        origin = Origin{Float64}()
        pbody = Sphere(radius_body1, mass_body1)
        cbody = Sphere(radius_body2, mass_body2)
        joint = JointConstraint(eval(joint_world_body1)(origin, pbody))
        bodies = [pbody, cbody]
        joints = [joint]

        # contact model
        collision = SphereSphereCollision{Float64,2,3,6}(
                szeros(3),
                szeros(3),
                pbody.shape.r,
                cbody.shape.r)

        if friction_type == :nonlinear
            friction_parameterization = SA{Float64}[
                    1.0  0.0
                    0.0  1.0
            ]   
            body_body = NonlinearContact{Float64,8}(friction_coefficient, friction_parameterization, collision)
        elseif friction_type == :linear 
            friction_parameterization = SA{Float64}[
                0.0  1.0
                0.0 -1.0
                1.0  0.0
            -1.0  0.0
            ]
            body_body = LinearContact{Float64,12}(friction_coefficient, friction_parameterization, collision)
        elseif friction_type == :impact
            friction_parameterization = szeros(Float64, 0, 2)
            collision = SphereSphereCollision{Float64,2,3,6}(
                        szeros(3),
                        szeros(3),
                        pbody.shape.r,
                        cbody.shape.r)
            body_body = ImpactContact{Float64,2}(friction_parameterization, collision)
        end

        contacts = [ContactConstraint((body_body, pbody.id, cbody.id), name=:body_body)]

        # z-drop (w/ gravity)
        mechanism = Mechanism(origin, bodies, joints, contacts,
                    gravity=gravity,
                    timestep=timestep)
        return mechanism
    end


    # z drop
    mech = get_two_body(; 
        gravity=-9.81)

    # initialize
    mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
    mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
    mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
    mech.bodies[2].state.v15 = [0.0, 0.0, 0.0]

    # initial distance
    d = distance(mech.contacts[1].model.collision,
        mech.bodies[1].state.x2, mech.bodies[1].state.q2,
        mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
    @test abs(d - 1.0) < 1.0e-6

    # body 1 contact point
    cp = contact_point(:parent, mech.contacts[1].model.collision,
        mech.bodies[1].state.x2, mech.bodies[1].state.q2,
        mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
    @test norm(cp - [0.0; 0.0; 0.5], Inf) < 1.0e-6

    # body 2 contact point
    cc = contact_point(:child, mech.contacts[1].model.collision,
        mech.bodies[1].state.x2, mech.bodies[1].state.q2,
        mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
    @test norm(cc - [0.0; 0.0; 1.5], Inf) < 1.0e-6

    # contact normal
    cn = contact_normal(mech.contacts[1].model.collision,
        mech.bodies[1].state.x2, mech.bodies[1].state.q2,
        mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
    @test norm(cn - [0.0 0.0 -1.0], Inf) < 1.0e-6

    # contact tangent
    ct = contact_tangent(mech.contacts[1].model.collision,
        mech.bodies[1].state.x2, mech.bodies[1].state.q2,
        mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
    @test norm(ct - [0.0 1.0 0.0; -1.0 0.0 0.0], Inf) < 1.0e-6

    # graph
    am = adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
    @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

    # simulate
    storage = simulate!(mech, 2.0,
        verbose=false,
        record=true)

    # test no interpentration
    @test norm(storage.x[2][end] - [0.0; 0.0; 1.0], Inf) < 1.0e-3

    # z velocity
    mech = get_two_body(; 
        gravity=0.0)

    # initialize
    mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
    mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
    mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
    mech.bodies[2].state.v15 = [0.0, 0.0, -5.0]

    # simulate
    storage = simulate!(mech, 2.0,
        verbose=false,
        record=true)

    # test no interpentration
    @test storage.x[2][end][end] > 1.0
end
