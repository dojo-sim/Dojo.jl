function get_pendulum(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    length=1,
    color=RGBA(1, 0, 0),
    springs=0,
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    spring_offset=szeros(1),
    orientation_offset=one(Quaternion),
    T=Float64)

    # mechanism
    origin = Origin{T}()
    
    body = Box(0.1, 0.1, length, mass; color, name=:pendulum)
    bodies = [body]

    joint = JointConstraint(Revolute(origin, body, X_AXIS;
        parent_vertex=1.1*Z_AXIS, child_vertex=0.5*Z_AXIS,
        rot_spring_offset=spring_offset, orientation_offset), 
        name=:joint)
    joints = [joint]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    zero_coordinates!(mechanism)

    # construction finished
    return mechanism
end

function initialize_pendulum!(mechanism::Mechanism;
    angle=0.7,
    angular_velocity=0)
    joint = mechanism.joints[1]
    set_minimal_coordinates_velocities!(mechanism, joint;
        xmin=[angle, angular_velocity])
end
