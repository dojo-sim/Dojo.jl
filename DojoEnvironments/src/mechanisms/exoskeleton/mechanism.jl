function get_exoskeleton(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:model,
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    limits=true,
    joint_limits=Dict([
        (:sAA, [0, 90]*π/180),
        (:sFE, [0, 90]*π/180),
        (:sIE, [-80, 25]*π/180),
        (:eFE, [-125,0]*π/180),]),
    keep_fixed_joints=false, 
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "dependencies/$(string(urdf)).urdf")
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
    initialize_exoskeleton!(mechanism)

    # construction finished
    return mechanism
end

function initialize_exoskeleton!(mechanism::Mechanism;
    joint_angles=[pi/2;pi/2-0.1;0;-0.1])

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    set_minimal_coordinates!(mechanism, get_joint(mechanism, :sAA), [joint_angles[1]])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :sFE), [joint_angles[2]])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :sIE), [joint_angles[3]])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :eFE), [joint_angles[4]])

    return
end
