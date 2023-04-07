function get_halfcheetah(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:halfcheetah,
    springs=0,
    dampers=0,
    parse_springs=true,
    parse_dampers=true,
    limits=true,
    joint_limits=Dict([
        (:bthigh, [-0.52,1.05]), 
        (:bshin, [-0.785,0.785]), 
        (:bfoot, [-0.400,0.785]), 
        (:fthigh, [-1,0.7]), 
        (:fshin, [-1.20,0.87]), 
        (:ffoot, [-0.5,0.5])]),
    keep_fixed_joints=true, 
    friction_coefficient=0.4,
    contact_feet=true,
    contact_body=true,
    T=Float64)
    
    # mechanism
    path = joinpath(@__DIR__, "../dependencies/$(String(urdf)).urdf")
    mechanism = Mechanism(path; floating=false, T,
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
        body_names = [:ffoot, :bfoot]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0;0; -0.5 * contact_bodies[1].shape.shapes[1].rh[2]],
            [0;0; -0.5 * contact_bodies[2].shape.shapes[1].rh[2]],
        ]
        contact_radii = [
            contact_bodies[1].shape.shapes[1].rh[1]
            contact_bodies[2].shape.shapes[1].rh[1]
        ]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    if contact_body
        body_names = getfield.(mechanism.bodies, :name)
        body_names = deleteat!(body_names, findall(x->(x==:ffoot||x==:bfoot||x==:torso),body_names))
        body_names = [:torso; :torso; :torso; body_names]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [[+0.5 * contact_bodies[1].shape.shapes[1].shapes[1].rh[2]; 0; 0]]
            [[-0.5 * contact_bodies[2].shape.shapes[1].shapes[1].rh[2]; 0; 0]]
            [[+0.5 * contact_bodies[3].shape.shapes[1].shapes[1].rh[2] + 0.214; 0; 0.1935]]
            [[0;0; -0.5 * contact_bodies[i].shape.shapes[1].rh[2]] for i = 4:n]
        ]
        contact_radii = [
            contact_bodies[1].shape.shapes[1].shapes[1].rh[1]
            contact_bodies[2].shape.shapes[1].shapes[1].rh[1]
            contact_bodies[3].shape.shapes[1].shapes[1].rh[1]
            [contact_bodies[i].shape.shapes[1].rh[1] for i = 4:n]
        ]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_halfcheetah!(mechanism)

    # construction finished
    return mechanism
end

function initialize_halfcheetah!(mechanism::Mechanism; 
    body_position=[0, 0], body_orientation=0)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint),
                 [body_position[1] + 0.576509, body_position[2], body_orientation + 0.02792])

    return
end
