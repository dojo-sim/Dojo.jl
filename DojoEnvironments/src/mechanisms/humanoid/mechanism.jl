function get_humanoid(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:humanoid, 
    springs=0, 
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.8, 
	contact_feet=true, 
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
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # contacts
    contacts = ContactConstraint{T}[]

    if contact_feet
        # feet contacts
        body_names = [:left_foot; :left_foot; :right_foot; :right_foot]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0.5 * contact_bodies[1].shape.shapes[1].shapes[1].rh[2]; 0; 0],
            [-0.5 * contact_bodies[2].shape.shapes[1].shapes[1].rh[2]; 0; 0],
            [0.5 * contact_bodies[3].shape.shapes[1].shapes[1].rh[2]; 0; 0],
            [-0.5 * contact_bodies[4].shape.shapes[1].shapes[1].rh[2]; 0; 0],
        ]
        contact_radii = [
            contact_bodies[1].shape.shapes[1].shapes[1].rh[1]
            contact_bodies[2].shape.shapes[1].shapes[1].rh[1]
            contact_bodies[3].shape.shapes[1].shapes[1].rh[1]
            contact_bodies[4].shape.shapes[1].shapes[1].rh[1]
        ]
        
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_humanoid!(mechanism)

    # construction finished
    return mechanism
end

function initialize_humanoid!(mechanism::Mechanism; 
    body_position=[0, 0, 1.33], body_orientation=one(Quaternion))

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_base), 
        [body_position; Dojo.rotation_vector(body_orientation)])

    return
end
