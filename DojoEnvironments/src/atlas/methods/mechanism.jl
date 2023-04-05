function get_atlas(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:atlas_simple,
    springs=0,
    dampers=0, 
    parse_springs=true, 
    parse_dampers=true, 
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=true, 
    friction_coefficient=0.8, 
    contact_feet=true, 
    contact_body=true, 
    T=Float64)

    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(string(urdf)).urdf")
    mechanism = Mechanism(path; floating=!(urdf==:armless), T,
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

    # contacts
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact_feet
        # Foot contact
        locations = [
            [-0.08; -0.04; 0.015],
            [+0.12; -0.02; 0.015],
            [-0.08; +0.04; 0.015],
            [+0.12; +0.02; 0.015],
            ]
        n = length(locations)
        normal = [Z_AXIS for i = 1:n]
        contact_radius = [0.025 for i = 1:n]
        friction_coefficients = friction_coefficient * ones(T, n)
        
        names = ["RR", "FR", "RL", "RR"]

        foot_contacts1 = contact_constraint(get_body(mechanism, :l_foot), normal; 
            friction_coefficient=friction_coefficients, 
            contact_origins=locations, 
            contact_radius, 
            names=[Symbol("l_" .* name) for name in names])
        foot_contacts2 = contact_constraint(get_body(mechanism, :r_foot), normal; 
            friction_coefficient=friction_coefficients, 
            contact_origins=locations, 
            contact_radius, 
            names=[Symbol("r_" .* name) for name in names])
        contacts = [contacts..., foot_contacts1..., foot_contacts2...]
    end

    if contact_body
        # Hands 
        location = [0; 0; 0]
        normal = Z_AXIS
        contact_radius = 0.06
        name = "hand"

        push!(contacts, contact_constraint(get_body(mechanism, :l_hand), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("l_" .* name)))
        push!(contacts, contact_constraint(get_body(mechanism, :r_hand), normal;
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("r_" .* name)))

        # knee 
        location = [0.025; 0; 0.175]
        normal = Z_AXIS
        contact_radius = 0.075
        name = "knee"

        push!(contacts, contact_constraint(get_body(mechanism, :l_lleg), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("l_" .* name)))
        push!(contacts, contact_constraint(get_body(mechanism, :r_lleg), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("r_" .* name)))

        # clavicals 
        location_l = [0; -0.05; -0.075]
        location_r = [0; -0.05; -0.075]

        normal = Z_AXIS
        contact_radius = 0.11
        name = "clav"

        push!(contacts, contact_constraint(get_body(mechanism, :l_clav), normal;
            friction_coefficient, 
            contact_origin=location_l, 
            contact_radius, 
            name=Symbol("l_" .* name)))
        push!(contacts, contact_constraint(get_body(mechanism, :r_clav), normal; 
            friction_coefficient, 
            contact_origin=location_r, 
            contact_radius, 
            name=Symbol("r_" .* name)))

        # pelvis
        location = [0; 0; 0.05]
        normal = Z_AXIS
        contact_radius = 0.19
        name = "pelvis"

        push!(contacts, contact_constraint(get_body(mechanism, :pelvis), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol(name)))

        # elbow 
        location = [0; -0.185; 0]
        normal = Z_AXIS
        contact_radius = 0.085
        name = "elbow"

        push!(contacts, contact_constraint(get_body(mechanism, :l_uarm), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("l_" .* name)))
        push!(contacts, contact_constraint(get_body(mechanism, :r_uarm), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol("r_" .* name)))

        # head
        location = [0; 0; 0]
        normal = Z_AXIS
        contact_radius = 0.175
        name = "head"

        push!(contacts, contact_constraint(get_body(mechanism, :head), normal;
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol(name)))

        # backpack
        location = [-0.095; 0; 0.25]
        normal = Z_AXIS
        contact_radius = 0.15
        name = "backpack_bottom"

        push!(contacts, contact_constraint(get_body(mechanism, :utorso), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol(name)))

        location = [-0.095; 0; -0.2]
        normal = Z_AXIS
        contact_radius = 0.15
        name = "backpack_top"

        push!(contacts, contact_constraint(get_body(mechanism, :utorso), normal; 
            friction_coefficient, 
            contact_origin=location, 
            contact_radius, 
            name=Symbol(name)))
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [0; 0; 0.9385; 0; 0; 0])

    # construction finished
    return mechanism
end

function initialize_atlas!(mechanism::Mechanism;
    model_type=:simple,
    body_position=[0, 0, 0.2],
    body_orientation=[0, 0, 0],
    link_linear_velocity=[zeros(3)  for i=1:length(mechanism.bodies)],
    link_angular_velocity=[zeros(3) for i=1:length(mechanism.bodies)],
    hip_orientation=0, 
    knee_orientation=0)

    body_position += (model_type == :armless) ? [0, 0, 0.9385 + 0.14853] : [0, 0, 0.9385]

    # positions
    try
        set_minimal_coordinates!(mechanism,
                get_joint(mechanism, :floating_base),
                [body_position; body_orientation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :back_bkxyz), [0, 0, 0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpxyz), [0, -hip_orientation, 0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpxyz), [0, -hip_orientation, 0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_kny), [knee_orienation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_kny), [knee_orienation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_akxy), [hip_orientation-knee_orienation, 0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_akxy), [hip_orientation-knee_orienation, 0])
    catch
        nothing
    end

    zero_velocity!(mechanism)

    return nothing
end

function initialize_atlas_stance!(mech::Mechanism;
    body_position=[0, 0, 0.2],
    body_orientation=[0, 0, 0],
    link_linear_velocity=[zeros(3)  for i=1:length(mechanism.bodies)],
    link_angular_velocity=[zeros(3) for i=1:length(mechanism.bodies)],
    hip_orientation=0, 
    knee_orienation=0)
    
    body_position += [0, 0, 0.9385]

    # positions
    try
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; body_orientation])

        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpxyz), [0, -hip_orientation, 0])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpxyz), [0, -hip_orientation, 0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_kny), [knee_orienation])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_kny), [knee_orienation])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_akxy), [hip_orientation - knee_orienation, 0])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_akxy), [hip_orientation - knee_orienation, 0])

        set_minimal_coordinates!(mechanism, get_joint(mechanism, :back_bkx), [0.0  * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :back_bky), [0.04 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :back_bkz), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_elx), [0.25 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_ely), [0.5 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_shx), [-0.5 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_shz), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_mwx), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_uwy), [0 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_arm_lwy), [0])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_akx), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_aky), [-0.1 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpx), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpy), [-0.1 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_hpz), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :l_leg_kny), [0.2 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :neck_ay), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_elx), [-0.25 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_ely), [0.5 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_shx), [0.5 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_shz), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_mwx), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_uwy), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_arm_lwy), [0 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_akx), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_aky), [-0.1 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpx), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpy), [-0.1 * π])
        # set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_hpz), [0 * π])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, :r_leg_kny), [0.2 * π])
    catch
        nothing
    end

    zero_velocity!(mechanism)

    return nothing
end