function get_rexhopper(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:rexhopper,
    springs=0.0,
    dampers=0.0,
    parse_springs=true, 
    parse_dampers=true,
    limits=true,
    joint_limits=Dict([
        (:joint1, [-0.7597,-1.8295]), 
        (:joint3, [1.8295,0.7597])])
    keep_fixed_joints=false, 
    friction_coefficient=1.0,
    contact_feet=true,
    contact_body=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(String(urdf)).urdf")
    mech = Mechanism(path; floating=true, T,
        gravity, timestep, input_scaling, 
        parse_damper, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mech.joints, springs)
    !parse_dampers && set_dampers!(mech.joints, dampers)

    # joint limits
    if limits
        joints = set_limits(mech, joint_limits)

        mech = Mechanism(Origin{T}(), mech.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = Origin{T}()
    bodies = mech.bodies
    joints = mech.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        normal = Z_AXIS

        link3 = get_body(mech, :link3)
        link2 = get_body(mech, :link2)
        foot_radius = 0.0203
        ankle_radius = 0.025
        base_radius = 0.14
        p = [0.1685; 0.0025; -0.0055]
        o = foot_radius
        push!(contacts, contact_constraint(link3, normal;
            friction_coefficient,
            contact_origin=p,
            contact_radius=o,
            contact_type,
            name=:foot))
        p = [-0.10; -0.002; 0.01]
        o = ankle_radius
        push!(contacts, contact_constraint(link3, normal;
            friction_coefficient,
            contact_origin=p,
            contact_radius=o,
            contact_type,
            name=:ankle3))
        p = [0.24; 0.007; 0.005]
        push!(contacts, contact_constraint(link2, normal;
            friction_coefficient,
            contact_origin=p,
            contact_radius=o,
            contact_type,
            name=:ankle2))
        base_link = get_body(mech, :base_link)
        pl = [0.0; +0.075; 0.03]
        pr = [0.0; -0.075; 0.03]
        o = base_radius
        push!(contacts, contact_constraint(base_link, normal;
            friction_coefficient,
            contact_origin=pl,
            contact_radius=o,
            contact_type,
            name=:torso_left))
        push!(contacts, contact_constraint(base_link, normal;
            friction_coefficient,
            contact_origin=pr,
            contact_radius=o,
            contact_type,
            name=:torso_right))

        # set_minimal_coordinates!(mech, get_joint(mech, :floating_base), [0,0,1.0, 0,0,0])
    end

    mech = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # construction finished
    return mech
end

function initialize_rexhopper!(mechanism::Mechanism{T};
        body_position=zeros(3),
        body_orientation=zeros(3),
        body_linear_velocity=zeros(3),
        body_angular_velocity=zeros(3)) where T

    zero_velocity!(mechanism)
    body_position += [0.0, 0.0, 0.3148]

    for joint in mechanism.joints
        !(joint.name in (:loop_joint, :floating_joint)) && set_minimal_coordinates!(mechanism, joint, zeros(input_dimension(joint)))
    end

    set_minimal_coordinates!(mechanism, get_joint(mechanism,
        :floating_base), [body_position; body_orientation])
    set_minimal_velocities!(mechanism, get_joint(mechanism,
        :floating_base), [body_linear_velocity; body_angular_velocity])
end
