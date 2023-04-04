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
    path = joinpath(@__DIR__, "../dependencies/$(String(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    # joint limits
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # zero configuration
    zero_coordinates!(mechanism)

    # construction finished
    return mechanism
end

function initialize_fourbar!(mechanism::Mechanism; 
    angle=0, 
    angular_velocity=szeros(2))

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :jointb1); 
        xmin=[ -angle, angular_velocity[1]])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :joint12); 
        xmin=[+2angle, 0])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :jointb3); 
        xmin=[ +angle, angular_velocity[2]])
    set_minimal_coordinates_velocities!(mechanism, get_joint(mechanism, :joint34); 
        xmin=[-2angle, 0])
end