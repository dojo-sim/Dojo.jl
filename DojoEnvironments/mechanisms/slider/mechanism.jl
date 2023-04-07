function get_slider(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    color=RGBA(1, 0, 0),
    springs=0, 
    dampers=0,
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    T=Float64)

    # mechanism
    origin = Origin{T}()

    pbody = Box(0.1, 0.1, 1, 1; color)
    bodies = [pbody]

    joint_between_origin_and_pbody = JointConstraint(Prismatic(origin, pbody, Z_AXIS; 
        child_vertex=Z_AXIS/2))
    joints = [joint_between_origin_and_pbody]

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

    # zero coordinates
    initialize_slider!(mechanism)

    # construction finished
    return mechanism
end

function initialize_slider!(mechanism::Mechanism; 
    position=mechanism.bodies[1].shape.xyz[3], velocity=0)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    set_minimal_coordinates!(mechanism, mechanism.joints[1], [position])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [velocity])

    return
end
