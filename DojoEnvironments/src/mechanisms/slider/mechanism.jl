function get_slider(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    color=RGBA(1, 0, 0),
    springs=0, 
    dampers=0,
    joint_limits=Dict(),
    keep_fixed_joints=true, 
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
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # zero coordinates
    initialize_slider!(mechanism)

    # construction finished
    return mechanism
end

function initialize_slider!(mechanism::Mechanism; 
    position=0, velocity=0)

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    child_vertex = mechanism.joints[1].translational.vertices[2]
    set_maximal_configurations!(mechanism.origin, mechanism.bodies[1];
        child_vertex=child_vertex - Z_AXIS*position)
    set_minimal_velocities!(mechanism, mechanism.joints[1], [velocity])

    return
end
