function get_quadruped(; 
    timestep=0.01, 
    gravity=[0.0; 0.0; -9.81], 
    friction_coefficient=0.8, 
    spring=0.0,
    damper=0.0, 
    contact_feet=true, 
    contact_body=true, 
    limits=true,
    spring_offset=true,
    path=joinpath(@__DIR__, "../deps/quadruped.urdf"),
    joint_limits=[[-0.5, -0.5, -2.5,],
                  [ 0.5,  1.5, -1.0,]],
    T=Float64)

    mech = Mechanism(path, true, T, 
        gravity=gravity, 
        timestep=timestep, 
        spring=spring, 
        damper=damper)

    # Adding springs and dampers
    for joint in mech.joints[2:end]
        joint.damper = true
        joint.spring = true
        joint.translational.spring=spring
        joint.translational.damper=damper
        joint.rotational.spring=spring
        joint.rotational.damper=damper
    end

    # joint limits
    joints = deepcopy(mech.joints)
    if limits
        for group in [:FR, :FL, :RR, :RL]
            hip_joint = get_joint(mech, Symbol(group, :_hip_joint))
            joints[hip_joint.id] = add_limits(mech, hip_joint,
                rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

            thigh_joint = get_joint(mech, Symbol(group, :_thigh_joint))
            joints[thigh_joint.id] = add_limits(mech, thigh_joint,
                rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

            calf_joint = get_joint(mech, Symbol(group, :_calf_joint))
            joints[calf_joint.id] = add_limits(mech, calf_joint,
                rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])
        end
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...], 
            gravity=gravity,
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end

    origin = Origin{T}()
    bodies = mech.bodies
    eqs = mech.joints

    contacts = ContactConstraint{T}[]
    if contact_feet
        # Foot contact
        normal = [0.0; 0.0; 1.0]
        foot_contact = [-0.006; 0.0; -0.092]
        foot_contact_radius = 0.021
        hip_contact = [0.0; 0.05; 0.0]
        hip_contact_radius = 0.05
        elbow_contactR = [-0.005; -0.023; -0.16]
        elbow_contactL = [-0.005; 0.023; -0.16]
        elbow_contact_radius = 0.023

        foot_contacts1 = contact_constraint(get_body(mech,:FR_calf), normal; 
            friction_coefficient=friction_coefficient, 
            contact_point=foot_contact, 
            contact_radius=foot_contact_radius, 
            name=:FR_contact)
        foot_contacts2 = contact_constraint(get_body(mech,:FL_calf), normal; 
            friction_coefficient=friction_coefficient, 
            contact_point=foot_contact, 
            contact_radius=foot_contact_radius, 
            name=:FL_contact)
        foot_contacts3 = contact_constraint(get_body(mech,:RR_calf), normal; 
            friction_coefficient=friction_coefficient, 
            contact_point=foot_contact, 
            contact_radius=foot_contact_radius, 
            name=:RR_contact)
        foot_contacts4 = contact_constraint(get_body(mech,:RL_calf), normal; 
            friction_coefficient=friction_coefficient, 
            contact_point=foot_contact, 
            contact_radius=foot_contact_radius, 
            name=:RL_contact)
        contacts = [contacts...; foot_contacts1; foot_contacts2; foot_contacts3; foot_contacts4]

        if contact_body
            elbow_contacts1 = contact_constraint(get_body(mech,:FR_thigh), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=elbow_contactR, 
                contact_radius=elbow_contact_radius, 
                name=:FR_hip_contact)
            elbow_contacts2 = contact_constraint(get_body(mech,:FL_thigh), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=elbow_contactL, 
                contact_radius=elbow_contact_radius, 
                name=:FL_hip_contact)
            elbow_contacts3 = contact_constraint(get_body(mech,:RR_thigh), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=elbow_contactR, 
                contact_radius=elbow_contact_radius, 
                name=:RR_hip_contact)
            elbow_contacts4 = contact_constraint(get_body(mech,:RL_thigh), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=elbow_contactL, 
                contact_radius=elbow_contact_radius, 
                name=:RL_hip_contact)
            push!(contacts, elbow_contacts1, elbow_contacts2, elbow_contacts3, elbow_contacts4)
            hip_contacts1 = contact_constraint(get_body(mech,:FR_hip), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=-hip_contact, 
                contact_radius=hip_contact_radius, 
                name=:FR_hip_contact)
            hip_contacts2 = contact_constraint(get_body(mech,:FL_hip), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=+hip_contact, 
                contact_radius=hip_contact_radius, 
                name=:FL_hip_contact)
            hip_contacts3 = contact_constraint(get_body(mech,:RR_hip), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=-hip_contact, 
                contact_radius=hip_contact_radius, 
                name=:RR_hip_contact)
            hip_contacts4 = contact_constraint(get_body(mech,:RL_hip), normal; 
                friction_coefficient=friction_coefficient, 
                contact_point=+hip_contact, 
                contact_radius=hip_contact_radius, 
                name=:RL_hip_contact)
            push!(contacts, hip_contacts1, hip_contacts2, hip_contacts3, hip_contacts4)
        end

        set_minimal_coordinates!(mech, get_joint(mech, :floating_base), [0.0; 0.0; 0.32; 0.0; 0.0; 0.0])
        
        mech = Mechanism(origin, bodies, eqs, contacts, 
            gravity=gravity, 
            timestep=timestep, 
            spring=spring, 
            damper=damper)
    end

    # Spring Offset
    if spring_offset
        θ_thigh = 0.9
        θ_calf = -1.425
        get_node(mech, :FR_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mech, :FL_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mech, :RR_hip_joint).rotational.spring_offset=0.0*sones(1)
        get_node(mech, :RL_hip_joint).rotational.spring_offset=0.0*sones(1)

        get_node(mech, :FR_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mech, :FL_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mech, :RR_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)
        get_node(mech, :RL_thigh_joint).rotational.spring_offset=θ_thigh*sones(1)

        get_node(mech, :FR_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mech, :FL_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mech, :RR_calf_joint).rotational.spring_offset=θ_calf*sones(1)
        get_node(mech, :RL_calf_joint).rotational.spring_offset=θ_calf*sones(1)
    end
    return mech
end

function initialize_quadruped!(mechanism::Mechanism{T}; 
    body_position=[0.0, 0.0, 0.0],
    body_orientation=[0.0, 0.0, 0.0], 
    body_linear_velocity=[0.0, 0.0, 0.0], 
    body_angular_velocity=[0.0, 0.0, 0.0], 
    link_angle=0.95) where T
    #TODO: add more leg angles

    body_position += [0.0, 0.0, 0.32]
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), 
        [body_position; body_orientation])

    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FR_hip_joint), [0.0])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FR_thigh_joint), [link_angle])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FR_calf_joint), [-1.5 * link_angle])

    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FL_hip_joint), [0.0])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FL_thigh_joint), [link_angle * 0.9])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :FL_calf_joint), [-1.5 * link_angle])

    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RR_hip_joint), [0.0])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RR_thigh_joint), [link_angle * 0.9])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RR_calf_joint), [-1.5 * link_angle])

    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RL_hip_joint), [0.0])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RL_thigh_joint), [link_angle])
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :RL_calf_joint), [-1.5 * link_angle])

    set_minimal_velocities!(mechanism, get_joint(mechanism, :floating_base), 
        [body_linear_velocity; body_angular_velocity])

    set_minimal_velocities!(mechanism, get_joint(mechanism, :FR_hip_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :FR_thigh_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :FR_calf_joint), [0.0])

    set_minimal_velocities!(mechanism, get_joint(mechanism, :FL_hip_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :FL_thigh_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :FL_calf_joint), [0.0])

    set_minimal_velocities!(mechanism, get_joint(mechanism, :RR_hip_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :RR_thigh_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :RR_calf_joint), [0.0])

    set_minimal_velocities!(mechanism, get_joint(mechanism, :RL_hip_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :RL_thigh_joint), [0.0])
    set_minimal_velocities!(mechanism, get_joint(mechanism, :RL_calf_joint), [0.0])
end