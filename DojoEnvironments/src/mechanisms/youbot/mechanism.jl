function get_youbot(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:youbot,
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict([
        (:arm_joint_1, [-2.95,2.95]), 
        (:arm_joint_2, [-1.57,1.13]), 
        (:arm_joint_3, [-2.55,2.55]), 
        (:arm_joint_4, [-1.78,1.78]), 
        (:arm_joint_5, [-2.92,2.92]), 
        # (:gripper_finger_joint_l, [0,0.03]), 
        # (:gripper_finger_joint_r, [-0.03,0]),
    ]),
    keep_fixed_joints=false,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "dependencies/$(String(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    # joint limits    
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_youbot!(mechanism)

    # construction finished
    return mechanism
end

function initialize_youbot!(mechanism;
    body_position=zeros(2), body_orientation=0, arm_angles=zeros(5))
    
    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    base_joint = get_joint(mechanism, :base_footprint_joint)
    arm_joint_1 = get_joint(mechanism, :arm_joint_1)
    arm_joint_2 = get_joint(mechanism, :arm_joint_2)
    arm_joint_3 = get_joint(mechanism, :arm_joint_3)
    arm_joint_4 = get_joint(mechanism, :arm_joint_4)
    arm_joint_5 = get_joint(mechanism, :arm_joint_5)

    base_joint !== nothing && set_minimal_coordinates!(mechanism, base_joint, [body_position;body_orientation])
    arm_joint_1 !== nothing && set_minimal_coordinates!(mechanism, get_joint(mechanism, :arm_joint_1), [arm_angles[1]])
    arm_joint_2 !== nothing && set_minimal_coordinates!(mechanism, get_joint(mechanism, :arm_joint_2), [arm_angles[2]])
    arm_joint_3 !== nothing && set_minimal_coordinates!(mechanism, get_joint(mechanism, :arm_joint_3), [arm_angles[3]])
    arm_joint_4 !== nothing && set_minimal_coordinates!(mechanism, get_joint(mechanism, :arm_joint_4), [arm_angles[4]])
    arm_joint_5 !== nothing && set_minimal_coordinates!(mechanism, get_joint(mechanism, :arm_joint_5), [arm_angles[5]])
    return
end