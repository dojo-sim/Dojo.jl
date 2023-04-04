function get_snake(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    num_bodies=2,
    length=1.0,
    radius=0.05,
    springs=0.0,
    dampers=0.0,
    limits=false,
    joint_limits=Dict(),
    joint_type=:Spherical,
    friction_coefficient=0.8,
    contact=true,
    contact_type=:nonlinear,   
    T=Float64)

    # mechanism
    origin = Origin{T}()

    bodies = [Box(3 * radius, 2 * radius, length, length, color=RGBA(0.9, 0.9, 0.9)) for i = 1:num_bodies]

    jointb1 = JointConstraint(Floating(origin, bodies[1]))

    joints = [
        jointb1;
        [
            JointConstraint(Prototype(joint_type, bodies[i - 1], bodies[i], X_AXIS;
            parent_vertex=-Z_AXIS*length/2, child_vertex=Z_AXIS*length/2)) for i = 2:num_bodies
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

        contacts1 = contact_constraint(bodies, normal;
            friction_coefficient,
            contact_origins=fill(Z_AXIS*length/2, n),
            contact_type) # we need to duplicate point for prismatic joint for instance
        contacts2 = contact_constraint(bodies, normal;
            friction_coefficient,
            contact_origins=fill(-Z_AXIS*length/2, n),
            contact_type)
        contacts = [contacts1; contacts2]        
    end

    mech = Mechanism(origin, bodies, joints, [contacts1; contacts2];
        gravity, timestep, input_scaling)

    # construction finished
    return mech
end

function initialize_snake!(mechanism::Mechanism{T};
    base_position=[0.0, -0.5, 0.0],
    base_orientation=RotX(0.6 * π),
    base_linear_velocity=zeros(3),
    base_angular_velocity=zeros(3),
    relative_linear_velocity=zeros(3),
    relative_angular_velocity=zeros(3)) where T

    pbody = mechanism.bodies[1]
    h = pbody.shape.xyz[3]
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
