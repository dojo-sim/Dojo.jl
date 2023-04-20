function get_uuv(; 
    timestep=0.01, 
    input_scaling=timestep,
    gravity=0, 
    urdf=:mini_tortuga_fixed_rotors,
    springs=0, 
    dampers=0, 
    parse_springs=true, 
    parse_dampers=true, 
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.5,
    contact_body=true,
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

    # contacts 
    contacts = ContactConstraint{T}[]

    if contact_body
        # base contact
        body_in_contact = get_body(mechanism, :base_link)
        n = 2
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_radii = fill(0.21,n)
        contact_origins = [
            [0.12; 0; 0.07], 
            [-0.12; 0; 0.07], 
        ]
        contacts = [contacts;contact_constraint(body_in_contact, normals; friction_coefficients, contact_radii, contact_origins)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_uuv!(mechanism)


    # construction finished
    return mechanism
end

function initialize_uuv!(mechanism::Mechanism; 
    body_position=Z_AXIS, body_orientation=one(Quaternion))
    
    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), [body_position; Dojo.rotation_vector(body_orientation)])
    
    return
end