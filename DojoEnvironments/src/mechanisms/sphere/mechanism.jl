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
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}(name=:origin)
    sphere  = Sphere(radius, mass; color, name=:sphere)
    bodies = [sphere]
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
    contacts = ContactConstraint{T}[]

    if contact
        contact_radius = radius
        contacts = [contacts;contact_constraint(sphere, Z_AXIS; friction_coefficient, contact_radius, contact_type)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_sphere!(mechanism)

    # construction finished
    return mechanism
end

function initialize_sphere!(mechanism::Mechanism;
    position=Z_AXIS/2, orientation=one(Quaternion),
    velocity=[1;0;0], angular_velocity=zeros(3))

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    position += Z_AXIS*mechanism.bodies[1].shape.r

    set_minimal_coordinates!(mechanism, mechanism.joints[1], [position; Dojo.rotation_vector(orientation)])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [velocity; angular_velocity])

    return
end
