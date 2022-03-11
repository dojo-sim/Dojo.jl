function get_snake(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,
    spring=0.0,
    damper=0.0,
    num_bodies=2,
    joint_type=:Spherical,
    length=1.0,
    radius=0.05,
    T=Float64)

    # Parameters
    ex = [1.0; 0.0; 0.0]
    ey = [0.0; 1.0; 0.0]
    ez = [0.0; 0.0; 1.0]

    vert11 = [0.0; 0.0; length / 2.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Box(3 * radius, 2 * radius, length, length, color=RGBA(0.9, 0.9, 0.9)) for i = 1:num_bodies]

    # Constraints
    jointb1 = JointConstraint(Floating(origin, bodies[1],
        spring=0.0,
        damper=0.0))

    if num_bodies > 1
        joints = [JointConstraint(Prototype(joint_type, bodies[i - 1], bodies[i], ex;
            parent_vertex=vert12,
            child_vertex=vert11,
            spring=spring,
            damper=damper)) for i = 2:num_bodies]
        joints = [jointb1; joints]
    else
        joints = [jointb1]
    end

    if contact
        n = num_bodies
        normal = [[0.0; 0.0; 1.0] for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts1 = contact_constraint(bodies, normal,
            friction_coefficient=friction_coefficient,
            contact_origins=fill(vert11, n),
            contact_type=contact_type) # we need to duplicate point for prismatic joint for instance
        contacts2 = contact_constraint(bodies, normal,
            friction_coefficient=friction_coefficient,
            contact_origins=fill(vert12, n),
            contact_type=contact_type)
        mech = Mechanism(origin, bodies, joints, [contacts1; contacts2],
            gravity=gravity,
            timestep=timestep,
            spring=spring,
            damper=damper)
    else
        mech = Mechanism(origin, bodies, joints,
            gravity=gravity,
            timestep=timestep,
            spring=spring,
            damper=damper)
    end
    return mech
end

function initialize_snake!(mechanism::Mechanism{T};
    base_position=[0.0, -0.5, 0.0],
    base_orientation=UnitQuaternion(RotX(0.6 * π)),
    base_linear_velocity=zeros(3),
    base_angular_velocity=zeros(3),
    relative_linear_velocity=zeros(3),
    relative_angular_velocity=zeros(3)) where T

    pbody = mechanism.bodies[1]
    h = pbody.shape.xyz[3]
    vert11 = [0.0; 0.0; h / 2.0]
    vert12 = -vert11

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody,
        child_vertex=base_position,
        Δq=base_orientation)
    set_maximal_velocities!(pbody,
        v=base_linear_velocity,
        ω=base_angular_velocity)

    previd = pbody.id
    for body in mechanism.bodies[2:end]
        set_maximal_configurations!(get_body(mechanism, previd), body,
            parent_vertex=vert12,
            child_vertex=vert11)
        set_maximal_velocities!(get_body(mechanism, previd), body,
            parent_vertex=vert12,
            child_vertex=vert11,
            Δv=relative_linear_velocity,
            Δω=relative_angular_velocity)
        previd = body.id
    end
end
