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
    contacts = ContactConstraint{T}[]

    if contact_feet
        # feet contacts
        left_names=[Symbol("l_" .* name) for name in ["RR", "FR", "RL", "RR"]]
        right_names=[Symbol("r_" .* name) for name in ["RR", "FR", "RL", "RR"]]
        n = length(left_names)
        left_foot = get_body(mechanism, :l_foot)
        right_foot = get_body(mechanism, :r_foot)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [-0.08; -0.04; 0.015],
            [+0.12; -0.02; 0.015],
            [-0.08; +0.04; 0.015],
            [+0.12; +0.02; 0.015],
        ]
        contact_radii = fill(0.025,n)
        
        contacts = [
            contacts
            contact_constraint(left_foot, normals; friction_coefficients, contact_origins, contact_radii, names=left_names)
            contact_constraint(right_foot, normals; friction_coefficients, contact_origins, contact_radii, names=right_names)
        ]
    end

    if contact_body
        body_names = [:l_hand, :r_hand, :l_lleg, :r_lleg, :l_clav, :r_clav, :pelvis, :l_uarm, :r_uarm, :head, :utorso, :utorso]
        names = [:l_hand, :r_hand, :l_knee, :r_knee, :l_clavis, :r_clavis, :pelvis, :l_elbow, :r_elbow, :head, :backpack_bottom, :backpack_top]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            zeros(3),           # l_hand
            zeros(3),           # r_hand
            [0.025; 0; 0.175],  # l_lleg
            [0.025; 0; 0.175],  # r_lleg
            [0; -0.05; -0.075], # l_clav
            [0; -0.05; -0.075], # r_clav
            [0; 0; 0.05],       # pelvis
            [0; -0.185; 0],     # l_uarm
            [0; -0.185; 0],     # r_uarm
            zeros(3),           # head
            [-0.095; 0; 0.25],  # backpack_bottom
            [-0.095; 0; -0.2],  # backpack_top
        ]
        contact_radii = [
            0.06,   # l_hand
            0.06,   # r_hand
            0.075,  # l_lleg
            0.075,  # r_lleg
            0.11,   # l_clav
            0.11,   # r_clav
            0.19,   # pelvis
            0.085,  # l_uarm
            0.085,  # r_uarm
            0.175,  # head
            0.15,   # backpack_bottom
            0.15,   # backpack_top
        ]
        contacts = [
            contacts
            contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii, names)
        ]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_atlas!(mechanism)

    # construction finished
    return mechanism
end

function initialize_atlas!(mechanism::Mechanism;
    body_position=[0, 0, 0.9385], body_orientation=one(Quaternion))

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; rotation_vector(body_orientation)])

    return
end