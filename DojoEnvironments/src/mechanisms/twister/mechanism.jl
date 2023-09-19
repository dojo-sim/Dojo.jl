function get_twister(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    num_bodies=5,
    height=1, 
    radius=0.05,
    color=RGBA(1, 0, 0),
    springs=0, 
    dampers=0, 
    joint_limits=Dict(),
    joint_type=:Prismatic, 
    keep_fixed_joints=true, 
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,  
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(height, 3 * radius, 2 * radius, height; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Floating(origin, bodies[1]))

    axes = [X_AXIS, Y_AXIS, Z_AXIS]
    joints = JointConstraint{T}[
        jointb1;
        [
            JointConstraint(Dojo.Prototype(joint_type, bodies[i - 1], bodies[i], axes[i%3+1]; 
            parent_vertex=-X_AXIS*height/2, child_vertex=X_AXIS*height/2)) for i = 2:num_bodies
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
        contact_bodies = [bodies[1];bodies]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [X_AXIS*height/2]  
            fill(-X_AXIS*height/2,n-1)
        ]
        contacts = [contacts;ContactConstraint(contact_type, contact_bodies, normals; friction_coefficients, contact_origins)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_twister!(mechanism)

    # construction finished
    return mechanism
end

function initialize_twister!(mechanism::Mechanism;
    base_position=zeros(3), base_orientation=one(Quaternion),
    base_linear_velocity=zeros(3), base_angular_velocity=zeros(3))

    zero_velocities!(mechanism)
    zero_coordinates!(mechanism)

    set_minimal_coordinates!(mechanism, mechanism.joints[1], [base_position; Dojo.rotation_vector(base_orientation)])
    set_minimal_velocities!(mechanism, mechanism.joints[1], [base_linear_velocity; base_angular_velocity])

    return
end
