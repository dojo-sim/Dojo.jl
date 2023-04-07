function get_fourbar(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:fourbar,
    springs=0, 
    dampers=0,
    parse_springs=true,
    parse_dampers=true,
    limits=false,
    joint_limits=Dict(),
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
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    initialize_fourbar!(mechanism)

    # construction finished
    return mechanism
end

function initialize_fourbar!(mechanism::Mechanism; 
    base_angle=0, inner_angle=pi/4)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)
    
    loop_joint_id = get_joint(mechanism, :joint24).id
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :jointb1), [base_angle+inner_angle]; exclude_ids = [loop_joint_id])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :jointb3), [base_angle-inner_angle]; exclude_ids = [loop_joint_id])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :joint12), [-2*inner_angle]; exclude_ids = [loop_joint_id])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :joint34), [2*inner_angle]; exclude_ids = [loop_joint_id])

    return
end