function get_twister(; 
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    num_bodies=5,
    height=1.0, 
    radius=0.05,
    springs=0.0, 
    dampers=0.0, 
    limits=false,
    joint_limits=Dict(),
    joint_type=:Prismatic, 
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,  
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(3 * radius, 2 * radius, height, height, color=RGBA(1.0, 0.0, 0.0)) for i = 1:num_bodies]

    jointb1 = JointConstraint(Floating(origin, bodies[1]))

    axes = [X_AXIS, Y_AXIS, Z_AXIS]
    joints = [
        jointb1;
        [
            JointConstraint(Prototype(joint_type, bodies[i - 1], bodies[i], axes[i%3+1]; 
            parent_vertex=-Z_AXIS*height/2, child_vertex=Z_AXIS*height/2)) for i = 2:num_bodies
        ]
    ]

    mech = Mechanism(origin, bodies, joints;
        gravity, timestep, input_scaling)

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
        n = num_bodies
        normal = [Z_AXIS for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)
        contacts1 = contact_constraint(bodies[1], normal[1]; 
        friction_coefficient=friction_coefficient[1], 
            contact_origin=Z_AXIS*height/2, 
            contact_type) # to avoid duplicating the contact points
        contacts2 = contact_constraint(bodies, normal;
            friction_coefficient, 
            contact_origins=fill(Z_AXIS*height/2, n), 
            contact_type)
        contacts = [contacts1; contacts2]
    end

    mech = Mechanism(origin, bodies, joints, contacts; 
        gravity, timestep, input_scaling)

    # construction finished
    return mech
end

function initialize_twister!(mechanism::Mechanism{T}; 
    base_position=[0.0, -1.0, 0.0],
    base_orientation=RotX(0.6 * π),
    base_linear_velocity=zeros(3), 
    base_angular_velocity=zeros(3),
    relative_linear_velocity=zeros(3), 
    relative_angular_velocity=zeros(3)) where T

    bodies = mechanism.bodies
    pbody = bodies[1]
    h = 1.0
    vert11 = [0.0; 0.0; h / 2.0]
    vert12 = -vert11

    # set position and velocities
    set_maximal_configurations!(mechanism.origin, pbody, 
        child_vertex=base_position, 
        Δq=base_orientation)
    set_maximal_velocities!(pbody, 
        v=base_linear_velocity, 
        ω=base_angular_velocity)

    previd = pbody.id
    for body in mechanism.bodies[2:end]
        set_maximal_configurations!(get_body(mechanism, previd), body, 
            parent_vertex=vert12, 
            child_vertex=vert11)
        set_maximal_velocities!(get_body(mechanism, previd), body, 
            parent_vertex=vert12, 
            child_vertex=vert11,
            Δv=relative_linear_velocity, 
            Δω=relative_angular_velocity)
        previd = body.id
    end
end
