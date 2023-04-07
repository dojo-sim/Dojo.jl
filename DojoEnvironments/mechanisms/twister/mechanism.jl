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
    limits=false,
    joint_limits=Dict(),
    joint_type=:Prismatic, 
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,  
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(height, 3 * radius, 2 * radius, height; color) for i = 1:num_bodies]

    jointb1 = JointConstraint(Floating(origin, bodies[1]))

    axes = [X_AXIS, Y_AXIS, Z_AXIS]
    joints = [
        jointb1;
        [
            JointConstraint(Prototype(joint_type, bodies[i - 1], bodies[i], axes[i%3+1]; 
            parent_vertex=-X_AXIS*height/2, child_vertex=X_AXIS*height/2)) for i = 2:num_bodies
        ]
    ]

    mechanism = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

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
        contact_bodies = [bodies[1];bodies]
        n = length(contact_bodies)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [X_AXIS*height/2]  
            fill(-X_AXIS*height/2,n-1)
        ]
        contacts = [contacts;contact_constraint(contact_bodies, normals; friction_coefficients, contact_origins, contact_type)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_twister!(mechanism)

    # construction finished
    return mechanism
end

function initialize_twister!(mechanism::Mechanism)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)

    return
end