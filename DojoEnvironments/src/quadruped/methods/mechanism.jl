function get_quadruped(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:gazebo_a1, 
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    spring_offset=true,
    limits=true,
    joint_limits=Dict(vcat([[
        (Symbol(group,:_hip_joint), [-0.5,0.5]), 
        (Symbol(group,:_thigh_joint), [-0.5,1.5]), 
        (Symbol(group,:_calf_joint), [-2.5,-1])] 
        for group in [:FR, :FL, :RR, :RL]]...)),
    keep_fixed_joints=true, 
    friction_coefficient=0.8,
    contact_feet=true,
    contact_body=true,
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mechanism = Mechanism(path; floating=true, T,
        gravity, timestep, input_scaling, 
        parse_dampers, keep_fixed_joints)

    # springs and dampers
    !parse_springs && set_springs!(mechanism.joints, springs)
    !parse_dampers && set_dampers!(mechanism.joints, dampers)

    if spring_offset
        θ_thigh = 0.9
        θ_calf = -1.425
        get_node(mechanism, :FR_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mechanism, :FL_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mechanism, :RR_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mechanism, :RL_hip_joint).rotational.spring_offset=0.0*sones(1)

        get_node(mechanism, :FR_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mechanism, :FL_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mechanism, :RR_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mechanism, :RL_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)

        get_node(mechanism, :FR_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mechanism, :FL_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mechanism, :RR_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mechanism, :RL_calf_joint).rotational.spring_offset=θ_calf*sones(1)
    end

    # joint limits
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        # Foot contact
        normal = Z_AXIS
        foot_contact = [-0.006; 0; -0.092]
        foot_contact_radius = 0.021

        foot_contacts1 = contact_constraint(get_body(mechanism,:FR_calf), normal;
            friction_coefficient,
            contact_origin=foot_contact,
            contact_radius=foot_contact_radius,
            name=:FR_contact)
        foot_contacts2 = contact_constraint(get_body(mechanism,:FL_calf), normal;
            friction_coefficient,
            contact_origin=foot_contact,
            contact_radius=foot_contact_radius,
            name=:FL_contact)
        foot_contacts3 = contact_constraint(get_body(mechanism,:RR_calf), normal;
            friction_coefficient,
            contact_origin=foot_contact,
            contact_radius=foot_contact_radius,
            name=:RR_contact)
        foot_contacts4 = contact_constraint(get_body(mechanism,:RL_calf), normal;
            friction_coefficient,
            contact_origin=foot_contact,
            contact_radius=foot_contact_radius,
            name=:RL_contact)
        contacts = [contacts...; foot_contacts1; foot_contacts2; foot_contacts3; foot_contacts4]
    end

    if contact_body
        normal = Z_AXIS
        hip_contact = [0; 0.05; 0]
        hip_contact_radius = 0.05
        elbow_contactR = [-0.005; -0.023; -0.16]
        elbow_contactL = [-0.005; 0.023; -0.16]
        elbow_contact_radius = 0.023

        elbow_contacts1 = contact_constraint(get_body(mechanism,:FR_thigh), normal;
            friction_coefficient,
            contact_origin=elbow_contactR,
            contact_radius=elbow_contact_radius,
            name=:FR_hip_contact)
        elbow_contacts2 = contact_constraint(get_body(mechanism,:FL_thigh), normal;
            friction_coefficient,
            contact_origin=elbow_contactL,
            contact_radius=elbow_contact_radius,
            name=:FL_hip_contact)
        elbow_contacts3 = contact_constraint(get_body(mechanism,:RR_thigh), normal;
            friction_coefficient,
            contact_origin=elbow_contactR,
            contact_radius=elbow_contact_radius,
            name=:RR_hip_contact)
        elbow_contacts4 = contact_constraint(get_body(mechanism,:RL_thigh), normal;
            friction_coefficient,
            contact_origin=elbow_contactL,
            contact_radius=elbow_contact_radius,
            name=:RL_hip_contact)
        push!(contacts, elbow_contacts1, elbow_contacts2, elbow_contacts3, elbow_contacts4)
        hip_contacts1 = contact_constraint(get_body(mechanism,:FR_hip), normal;
            friction_coefficient,
            contact_origin=-hip_contact,
            contact_radius=hip_contact_radius,
            name=:FR_hip_contact)
        hip_contacts2 = contact_constraint(get_body(mechanism,:FL_hip), normal;
            friction_coefficient,
            contact_origin=+hip_contact,
            contact_radius=hip_contact_radius,
            name=:FL_hip_contact)
        hip_contacts3 = contact_constraint(get_body(mechanism,:RR_hip), normal;
            friction_coefficient,
            contact_origin=-hip_contact,
            contact_radius=hip_contact_radius,
            name=:RR_hip_contact)
        hip_contacts4 = contact_constraint(get_body(mechanism,:RL_hip), normal;
            friction_coefficient,
            contact_origin=+hip_contact,
            contact_radius=hip_contact_radius,
            name=:RL_hip_contact)
        push!(contacts, hip_contacts1, hip_contacts2, hip_contacts3, hip_contacts4)
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_quadruped!(mechanism)

    # construction finished
    return mechanism
end

function initialize_quadruped!(mechanism::Mechanism;
    body_position=[0, 0, 0], body_orientation=one(Quaternion),
    hip_angle=0, thigh_angle=pi/4, calf_angle=-pi/2)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    body_position += [0, 0, 0.43]
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; rotation_vector(body_orientation)])
    for group in [:FR, :FL, :RR, :RL]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol(group, :_hip_joint)), [hip_angle])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol(group, :_thigh_joint)), [thigh_angle])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol(group, :_calf_joint)), [calf_angle])
    end

    return
end
