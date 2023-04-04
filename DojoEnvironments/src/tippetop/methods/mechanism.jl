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
    friction_coefficient=0.4,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}(name=:origin)
    
    bodies = [
        Sphere(radius, mass; name=:sphere1, color=gray_light),
        Sphere(radius*scale, mass*scale^3; name=:sphere2, color=gray_light)]
    bodies[1].inertia = Diagonal([1.9, 2.1, 2])

    joints = [
        JointConstraint(Floating(origin, bodies[1]);
                name=:floating_joint),
        JointConstraint(Fixed(bodies[1], bodies[2];
                parent_vertex=[0,0,radius]),
                name = :fixed_joint)]

    mechanism = Mechanism(origin, bodies, joints;
        timestep, gravity, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

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

    if contact
        contact_origin = [0, 0, 0]
        normal = Z_AXIS
        contacts = [
            contact_constraint(get_body(mechanism, :sphere1), normal;
                friction_coefficient,
                contact_origin, 
                contact_radius=radius,
                contact_type),
            contact_constraint(get_body(mechanism, :sphere2), normal;
                friction_coefficient,
                contact_origin,
                contact_radius=radius * scale,
                contact_type)
            ]
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, get_joint(mechanism, :floating_joint), [0; 0; radius; 0; 0; 0])

    # construction finished
    return mechanism
end

function initialize_tippetop!(mechanism::Mechanism{T};
    body_position=zeros(3),
    body_orientation=zeros(3),
    body_linear_velocity=zeros(3),
    body_angular_velocity=zeros(3)) where T

    floating_joint = get_joint(mechanism, :floating_joint)
    fixed_joint = get_joint(mechanism, :fixed_joint)
    radius = fixed_joint.translational.vertices[1][3]

    zero_velocity!(mechanism)
    set_minimal_coordinates_velocities!(mechanism, floating_joint,
        xmin=[body_position; body_orientation; body_linear_velocity; body_angular_velocity])
    set_minimal_coordinates_velocities!(mechanism, fixed_joint)
    return nothing
end
