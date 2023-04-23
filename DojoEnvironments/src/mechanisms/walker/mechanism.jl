function get_walker(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:walker, 
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict([
        (:thigh, [0,150] * π / 180), 
        (:leg, [0,150] * π / 180), 
        (:foot, [-45,45] * π / 180), 
        (:thigh_left, [0,150] * π / 180), 
        (:leg_left, [0,150] * π / 180), 
        (:foot_left, [-45,45] * π / 180)]),
    keep_fixed_joints=false, 
    friction_coefficient=0.5,
    contact_feet=true,
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

    if contact_feet
        body_names = [:foot; :foot; :foot_left; :foot_left]
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0, 0, +0.5 * contact_bodies[1].shape.shapes[1].rh[2]],
            [0, 0, -0.5 * contact_bodies[2].shape.shapes[1].rh[2]],
            [0, 0, +0.5 * contact_bodies[3].shape.shapes[1].rh[2]],
            [0, 0, -0.5 * contact_bodies[4].shape.shapes[1].rh[2]],
        ]
        contact_radii = [
            contact_bodies[1].shape.shapes[1].rh[1]
            contact_bodies[2].shape.shapes[1].rh[1]
            contact_bodies[3].shape.shapes[1].rh[1]
            contact_bodies[4].shape.shapes[1].rh[1]
        ]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end
    
    if contact_body
        body_names = getfield.(mechanism.bodies, :name)
        body_names = deleteat!(body_names, findall(x->(x==:foot||x==:foot_left),body_names))
        contact_bodies = [get_body(mechanism, name) for name in body_names]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [[0;0; 0.5 * contact_bodies[i].shape.shapes[1].rh[2]] for i = 1:n]
        contact_radii = [contact_bodies[i].shape.shapes[1].rh[1] for i = 1:n]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_radii)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts; 
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_walker!(mechanism)

    # construction finished
    return mechanism
end

function initialize_walker!(mechanism::Mechanism; 
    body_position=[0, 0], body_orientation=0)

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)
    
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint),
        [body_position[1] + 1.25 , body_position[2], body_orientation])

    return
end
