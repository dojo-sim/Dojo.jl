function get_pendulum(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    link_length=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0,
    joint_limits=Dict(),
    spring_offset=szeros(1),
    orientation_offset=one(Quaternion),
    T=Float64)

    # mechanism
    origin = Origin{T}()
    
    body = Box(0.1, 0.1, link_length, mass; color, name=:pendulum)
    bodies = [body]

    joint = JointConstraint(Revolute(origin, body, X_AXIS;
        parent_vertex=(link_length+0.1)*Z_AXIS, child_vertex=0.5*link_length*Z_AXIS,
        rot_spring_offset=spring_offset, orientation_offset), 
        name=:joint)
    joints = [joint]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_pendulum!(mechanism)

    # construction finished
    return mechanism
end

function initialize_pendulum!(mechanism::Mechanism;
    angle=pi/4, angular_velocity=0)

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, mechanism.joints[1], [angle])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [angular_velocity])

    return
end
