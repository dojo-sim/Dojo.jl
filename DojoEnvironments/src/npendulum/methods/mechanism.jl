function get_npendulum(;
    timestep=0.01,
    input_scaling=timestep,
    gravity=-9.81,
    num_bodies=5,
    mass=1,
    length=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    base_joint_type=:Revolute,
    rest_joint_type=:Revolute,
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(0.05, 0.05, length, mass; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Prototype(base_joint_type, origin, bodies[1], X_AXIS;
        parent_vertex=Z_AXIS*num_bodies, child_vertex=Z_AXIS*length/2))

    joints = [
        jointb1;
        [
            JointConstraint(Prototype(rest_joint_type, bodies[i - 1], bodies[i], X_AXIS;
            parent_vertex=-Z_AXIS*length/2, child_vertex=Z_AXIS*length/2)) for i = 2:num_bodies
        ]
    ]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    zero_coordinates!(mechanism)

    # construction finished
    return mechanism
end

function initialize_npendulum!(mechanism::Mechanism{T};
    base_angle=π / 4,
    base_angular_velocity=[0, 0, 0],
    relative_linear_velocity=[0, 0, 0],
    relative_angular_velocity=[0, 0, 0]) where T

    pbody = mechanism.bodies[1]
    joint = mechanism.joints[1]
    vert11 = joint.translational.vertices[2]
    vert12 = -vert11

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody,
        child_vertex=vert11,
        Δq=RotX(base_angle))
    set_maximal_velocities!(pbody,
        ω=base_angular_velocity)

    previd = pbody.id
    for (i, body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        set_maximal_configurations!(get_body(mechanism, previd), body,
            parent_vertex=vert12,
            child_vertex=vert11)
        set_maximal_velocities!(get_body(mechanism, previd), body,
            parent_vertex=vert12,
            child_vertex=vert11,
            Δv=relative_linear_velocity, Δω=1 / i * relative_angular_velocity)
        previd = body.id
    end
end
