function get_youbot(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:youbot,
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
    initialize_youbot!(mechanism)

    # construction finished
    return mechanism
end

function initialize_youbot!(mechanism)
    
    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    return
end