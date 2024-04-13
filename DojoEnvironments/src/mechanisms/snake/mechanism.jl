function get_snake(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    num_bodies=2,
    link_length=1,
    radius=0.05,
    color=RGBA(0.9, 0.9, 0.9),
    springs=0,
    dampers=0,
    joint_limits=Dict(),
    joint_type=:Spherical,
    keep_fixed_joints=true, 
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,   
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(link_length, 3 * radius, 2 * radius, link_length; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Floating(origin, bodies[1]))

    joints = JointConstraint{T}[
        jointb1;
        [
            JointConstraint(Dojo.Prototype(joint_type, bodies[i - 1], bodies[i], X_AXIS;
            parent_vertex=-X_AXIS*link_length/2, child_vertex=X_AXIS*link_length/2)) for i = 2:num_bodies
        ]
    ]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

    # springs and dampers
    set_springs!(mechanism.joints, springs)
    set_dampers!(mechanism.joints, dampers)

    # joint limits    
    joints = set_limits(mechanism, joint_limits)
    mechanism = Mechanism(mechanism.origin, mechanism.bodies, joints;
        gravity, timestep, input_scaling)

    # contacts
    contacts = ContactConstraint{T}[]

    if contact
        contact_bodies = [bodies;bodies] # we need to duplicate contacts for prismatic joint for instance
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            fill(X_AXIS*link_length/2, Int64(n/2))
            fill(-X_AXIS*link_length/2, Int64(n/2))
        ]
        contacts = [contacts;ContactConstraint(contact_type, contact_bodies, normals; friction_coefficients, contact_origins)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_snake!(mechanism)

    # construction finished
    return mechanism
end

function initialize_snake!(mechanism::Mechanism;
    base_position=zeros(3), base_orientation=one(Quaternion),
    base_linear_velocity=zeros(3), base_angular_velocity=zeros(3))

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    set_minimal_coordinates!(mechanism, mechanism.joints[1], [base_position; Dojo.rotation_vector(base_orientation)])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [base_linear_velocity; base_angular_velocity])

    return
end
