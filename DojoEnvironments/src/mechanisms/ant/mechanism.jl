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
    path = joinpath(@__DIR__, "dependencies/$(string(urdf)).urdf")
    mechanism = Mechanism(path; floating=true, T,
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
        body_names = [:front_left_foot, :front_right_foot, :left_back_foot, :right_back_foot]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0.2; 0.2; 0], 
            [-0.2; 0.2; 0], 
            [-0.2; -0.2; 0], 
            [0.2; -0.2; 0]
        ]
        contact_radii = [body.shape.shapes[1].rh[1] for body in contact_bodies]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    if contact_body
        # torso contact
        torso = get_body(mechanism, :torso)
        contact_radius = torso.shape.r
        contacts = [contacts;contact_constraint(torso, Z_AXIS; friction_coefficient, contact_radius)]

        # elbow contact
        body_names = [:aux_1, :aux_2, :aux_3, :aux_4]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            -[0.1; 0.1; 0], 
            -[-0.1; 0.1; 0], 
            -[-0.1;-0.1; 0], 
            -[0.1; -0.1; 0]
        ]
        contact_radii = [body.shape.shapes[1].rh[1] for body in contact_bodies]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_ant!(mechanism)


    # construction finished
    return mechanism
end

function initialize_ant!(mechanism::Mechanism; 
    body_position=0.5*Z_AXIS, body_orientation=one(Quaternion), ankle_angle=0.25)
    
    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; Dojo.rotation_vector(body_orientation)])

    for i in [1, 4]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("hip_$i")), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("ankle_$i")), [ankle_angle * π])
    end

    for i in [2, 3]
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("hip_$i")), [0])
        set_minimal_coordinates!(mechanism, get_joint(mechanism, Symbol("ankle_$i")), [-ankle_angle * π])
    end
    
    return
end