function get_sphere(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    radius=0.5,
    color=RGBA(0.9, 0.9, 0.9),
    springs=0,
    dampers=0, 
    limits=false,
    joint_limits=Dict(),
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}(name=:origin)
    bodies = [Sphere(radius, mass; color, name=:sphere)]
    joints = [JointConstraint(Floating(origin, bodies[1]); name=:floating_joint)]

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
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact
        contact = [0, 0, 0]
        normal = Z_AXIS
        contacts = [contact_constraint(get_body(mechanism, :sphere), normal;
            friction_coefficient,
            contact_origin=contact,
            contact_radius=radius,
            contact_type)]
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, joints[1], [0,0,radius, 0,0,0])

    # construction finished
    return mechanism
end

function initialize_sphere!(mechanism::Mechanism{T};
    position=zeros(3),
    orientation=one(Quaternion),
    velocity=zeros(3),
    angular_velocity=zeros(3)) where T

    r = mechanism.bodies[1].shape.r
    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [position + [0, 0, r]; rotation_vector(orientation)])
    set_minimal_velocities!(mechanism, joint, [velocity; angular_velocity])
end
