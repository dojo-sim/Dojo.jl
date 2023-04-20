function get_panda(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:panda_end_effector,
    springs=0,
    dampers=5,
    parse_springs=true, 
    parse_dampers=false,
    joint_limits=Dict([
        (:joint1, [-2.8973, 2.8973]),
        (:joint2, [-1.7628, 1.7628]),
        (:joint3, [-2.8973, 2.8973]),
        (:joint4, [-3.0718,-0.0698]),
        (:joint5, [-2.8973, 2.8973]),
        (:joint6, [-0.0175, 3.7525]),
        (:joint7, [-2.8973, 2.8973]),
        (:jointf1, [-0.00, 0.04]),
        (:jointf2, [-0.00, 0.04])]),
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
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_panda!(mechanism)

    # construction finished
    return mechanism
end

function initialize_panda!(mechanism::Mechanism;
    joint_angles=[0;0.5;0;-0.5;0;0.5;0; zeros(input_dimension(mechanism)-7)],
    joint_velocities=zeros(input_dimension(mechanism)))

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    nu = input_dimension(mechanism)
    x = vcat([[joint_angles[i], joint_velocities[i]] for i=1:nu]...)
    set_minimal_state!(mechanism, x)

    return
end
