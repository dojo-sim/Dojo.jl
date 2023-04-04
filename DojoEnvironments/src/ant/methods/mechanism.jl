function get_ant(; 
    timestep=0.05, 
    input_scaling=timestep,
    gravity=-9.81, 
    urdf=:ant,
    springs=0, 
    dampers=0, 
    parse_springs=true, 
    parse_dampers=true, 
    limits=true,
    joint_limits=Dict([
        (:hip_1, [-30,30] * π / 180), 
        (:ankle_1, [30,70] * π / 180), 
        (:hip_2, [-30,30] * π / 180), 
        (:ankle_2, [-70,-30] * π / 180), 
        (:hip_3, [-30,30] * π / 180), 
        (:ankle_3, [-70,-30] * π / 180), 
        (:hip_4, [-30,30] * π / 180), 
        (:ankle_4, [30,70] * π / 180)]),
    keep_fixed_joints=true, 
    friction_coefficient=0.5,
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

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(Origin{T}(), mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts 
    origin = Origin{T}()
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        # foot contact
        normal = Z_AXIS
        foot_names = [:front_left_foot, :front_right_foot, :left_back_foot, :right_back_foot]
        foot = [get_body(mechanism, name) for name in foot_names]
        p = [[0.2; 0.2; 0], [-0.2; 0.2; 0], [-0.2; -0.2; 0], [0.2; -0.2; 0]]
        o = [f.shape.shapes[1].rh[1] for f in foot]
        contacts = [contact_constraint(foot[i], normal; friction_coefficient, contact_origin=p[i], contact_radius=o[i]) for i = 1:length(foot_names)]
    end

    if contact_body
        # torso contact
        torso = get_body(mechanism, :torso)
        p = [0; 0; 0]
        o = torso.shape.r
        torso_contacts = contact_constraint(torso, normal; friction_coefficient, contact_origin=p, contact_radius=o)

        # elbow contact
        elbow_names = [:aux_1, :aux_2, :aux_3, :aux_4]
        elbow = [get_body(mechanism, e) for e in elbow_names]
        p = [-[0.1; 0.1; 0], -[-0.1; 0.1; 0], -[-0.1;-0.1; 0], -[0.1; -0.1; 0]]
        o = [e.shape.shapes[1].rh[1] for e in elbow]
        elbow_contacts = [contact_constraint(elbow[i], normal; friction_coefficient, contact_origin=p[i], contact_radius=o[i]) for i = 1:length(elbow_names)]

        contacts = [contacts..., torso_contacts, elbow_contacts...]
    end

    mechanism = Mechanism(origin, bodies, joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [0; 0; 0.4; 0; 0; 0])


    # construction finished
    return mechanism
end

function initialize_ant!(mechanism::Mechanism; 
    body_position=[0; 0; 0.63], 
    body_orientation=[0; 0; 0 * π],
    ankle_orientation=0.25)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; body_orientation])

    for i in [1, 4]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("hip_$i")), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("ankle_$i")), [ankle_orientation * π])
    end

    for i in [2, 3]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("hip_$i")), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("ankle_$i")), [-ankle_orientation * π])
    end

    zero_velocity!(mechanism)
end