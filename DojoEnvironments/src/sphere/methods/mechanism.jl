function get_sphere(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1
    radius=0.5,
    springs=0.0,
    dampers=0.0, 
    limits=false,
    joint_limits=Dict(),
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}(name=:origin)
    bodies = [Sphere(radius, mass, name=:sphere)]
    joints = [JointConstraint(Floating(origin, bodies[1]), name=:floating_joint)]

    mech = Mechanism(origin, bodies, joints;
        timestep, gravity, input_scaling)

    # springs and dampers
    set_springs!(mech.joints, springs)
    set_dampers!(mech.joints, dampers)

    # joint limits    
    if limits
        joints = set_limits(mech, joint_limits)

        mech = Mechanism(Origin{T}(), mech.bodies, joints;
            gravity, timestep, input_scaling)
    end

    # contacts
    origin = Origin{T}()
    bodies = mech.bodies
    joints = mech.joints
    contacts = ContactConstraint{T}[]

    if contact
        contact = [0.0, 0.0, 0.0]
        normal = Z_AXIS
        contacts = [contact_constraint(get_body(mech, :sphere), normal;
            friction_coefficient,
            contact_origin=contact,
            contact_radius=radius,
            contact_type)]
        # set_minimal_coordinates!(mech, get_joint(mech, :floating_joint), [0.0; 0.0; radius; zeros(3)])
    end

    mech = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # construction finished
    return mech
end

function initialize_sphere!(mechanism::Mechanism{T};
    position=zeros(3),
    orientation=one(Quaternion),
    velocity=zeros(3),
    angular_velocity=zeros(3)) where T

    r = mechanism.bodies[1].shape.r
    joint = get_joint(mechanism, :floating_joint)
    zero_velocity!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [position + [0.0, 0.0, r]; rotation_vector(orientation)])
    set_minimal_velocities!(mechanism, joint, [velocity; angular_velocity])
end
