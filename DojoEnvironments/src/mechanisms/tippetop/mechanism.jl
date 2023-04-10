function get_tippetop(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    radius=0.5,
    scale=0.2,
    springs=0,
    dampers=0, 
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.4,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}(name=:origin)
    
    bodies = [
        Sphere(radius, mass; name=:sphere1, color=gray_light),
        Sphere(radius*scale, mass*scale^3; name=:sphere2, color=gray_light)
    ]
    bodies[1].inertia = Diagonal([1.9, 2.1, 2])

    joints = [
        JointConstraint(Floating(origin, bodies[1]); name=:floating_joint),
        JointConstraint(Fixed(bodies[1], bodies[2];
            parent_vertex=[0,0,radius]), name = :fixed_joint)
    ]

    mechanism = Mechanism(origin, bodies, joints;
        timestep, gravity, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mechanism, joint_limits)

        mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
            gravity, timestep, input_scaling)
    end
    
    # contacts
    contacts = ContactConstraint{T}[]

    if contact
        n = length(bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_radii = [radius;radius*scale]
        contacts = [contacts;contact_constraint(bodies, normals; friction_coefficients, contact_radii, contact_type)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_tippetop!(mechanism)

    # construction finished
    return mechanism
end

function initialize_tippetop!(mechanism::Mechanism{T};
    body_position=2*Z_AXIS*mechanism.bodies[1].shape.r, body_orientation=one(Quaternion),
    body_linear_velocity=zeros(3), body_angular_velocity=[0.0, 0.1, 50.0]) where T

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    floating_joint = mechanism.joints[1]
    set_minimal_coordinates!(mechanism, floating_joint, [body_position; rotation_vector(body_orientation)])
    set_minimal_velocities!(mechanism, floating_joint, [body_linear_velocity; body_angular_velocity])
end
