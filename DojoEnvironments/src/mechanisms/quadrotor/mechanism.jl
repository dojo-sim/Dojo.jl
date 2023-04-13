function get_quadrotor(; 
    timestep=0.01, 
    input_scaling=timestep,
    gravity=-9.81, 
    urdf=:pelican_fixed_rotors,
    springs=0, 
    dampers=0, 
    parse_springs=true, 
    parse_dampers=true, 
    limits=true,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.5,
    contact_rotors=true, 
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

    if contact_rotors
        # rotor contacts
        body_in_contact = get_body(mechanism, :base_link)
        n = 4
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_radii = fill(0.07,4) # smaller than rotor size to allow for pitching close to ground
        contact_origins = [
            [0.21; 0; 0.045], 
            [-0.21; 0; 0.045], 
            [0; 0.21; 0.045], 
            [0; -0.21; 0.045]
        ]
        contacts = [contacts;contact_constraint(body_in_contact, normals; friction_coefficients, contact_radii, contact_origins)]
    end

    if contact_body
        # base contact
        body_in_contact = get_body(mechanism, :base_link)
        n = 4
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0.11; 0; -0.085], 
            [-0.11; 0; -0.085], 
            [0; 0.11; -0.085], 
            [0; -0.11; -0.085]
        ]
        contacts = [contacts;contact_constraint(body_in_contact, normals; friction_coefficients, contact_origins)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_quadrotor!(mechanism)


    # construction finished
    return mechanism
end

function initialize_quadrotor!(mechanism::Mechanism; 
    body_position=0.085*Z_AXIS, body_orientation=one(Quaternion))
    
    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; Dojo.rotation_vector(body_orientation)])
    
    return
end