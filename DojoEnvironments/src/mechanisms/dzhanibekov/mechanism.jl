function get_dzhanibekov(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=0,
    springs=0,
    dampers=0, 
    color=RGBA(0.9,0.9,0.9,1),
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    T=Float64)

    # mechanism
    origin = Origin{Float64}()
    main_body = Capsule(0.1, 1, 1; color, name=:main)
    main_body.inertia = Diagonal([3e-2, 1e-3, 1e-1])
    side_body = Capsule(0.5 * 0.1, 0.35 * 1, 0.5 * 1;
        orientation_offset=RotY(0.5 * π),
        color=color, name=:side)
    bodies = [main_body, side_body]

    joint_float = JointConstraint(Floating(origin, bodies[1]); name=:floating)
    joint_fixed = JointConstraint(Fixed(bodies[1], bodies[2];
        parent_vertex=szeros(3), child_vertex=[-0.25 * 1; 0; 0]), name=:fixed)
    joints = [joint_float, joint_fixed]
    
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
    initialize_dzhanibekov!(mechanism)

    # construction finished
    return mechanism
end

function initialize_dzhanibekov!(mechanism::Mechanism;
    linear_velocity=zeros(3), angular_velocity=[10.0; 0.01; 0.0])

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)
    
    joint = mechanism.joints[1]
    set_minimal_coordinates!(mechanism, joint, [Z_AXIS; zeros(3)])
    set_minimal_velocities!(mechanism, joint, [linear_velocity; angular_velocity])

    return
end
