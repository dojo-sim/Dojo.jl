function get_pendulum(;
    timestep=0.01,
    gravity=-9.81,
    mass=1.0,
    len=1.0,
    spring=0.0,
    damper=0.0,
    spring_offset=szeros(1),
    axis_offset=one(Quaternion),
    T=Float64)

    # Parameters
    joint_axis = [1.0; 0.0; 0.0]
    width, depth = 0.1, 0.1
    parent_vertex = [0.0; 0.0; 1.1 * len] # joint connection point
    child_vertex = [0.0; 0.0; 0.5 * len] # joint connection point

    # Links
    origin = Origin{T}()
    body = Box(width, depth, len, mass, name=:pendulum)

    # Constraints
    joint = JointConstraint(Revolute(origin, body, joint_axis;
        parent_vertex=parent_vertex,
        child_vertex=child_vertex,
        spring=spring,
        damper=damper,
        rot_spring_offset=spring_offset,
        axis_offset=axis_offset), name=:joint)

    bodies = [body]
    joints = [joint]

    mech = Mechanism(origin, bodies, joints,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    return mech
end

function initialize_pendulum!(mechanism::Mechanism;
    angle=0.7,
    angular_velocity=0.0) where T
    joint = mechanism.joints[1]
    set_minimal_coordinates_velocities!(mechanism, joint;
        xmin=[angle, angular_velocity])
end

function get_npendulum(;
    timestep=0.01,
    gravity=-9.81,
    mass=1.0,
    len=1.0,
    spring=0.0,
    damper=0.0,
    num_bodies=5,
    basetype=:Revolute,
    joint_type=:Revolute,
    T=Float64)

    # Parameters
    ex = [1.0; 0.0; 0.0]
    r = 0.05
    vert11 = [0.0; 0.0; len / 2.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Box(r, r, len, mass, color=RGBA(1.0, 0.0, 0.0)) for i = 1:num_bodies]

    # Constraints
    jointb1 = JointConstraint(Prototype(basetype, origin, bodies[1], ex;
        child_vertex=vert11,
        spring=spring,
        damper=damper))
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

    mech = Mechanism(origin, bodies, joints,
        gravity=gravity,
        timestep=timestep)
    return mech
end

function initialize_npendulum!(mechanism::Mechanism{T};
    base_angle=π / 4.0,
    base_angular_velocity=[0.0, 0.0, 0.0],
    relative_linear_velocity=[0.0, 0.0, 0.0],
    relative_angular_velocity=[0.0, 0.0, 0.0]) where T

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
