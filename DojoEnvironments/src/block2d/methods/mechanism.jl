function get_block2d(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    edge_length=0.5,
    color=RGBA(1, 1, 0.),
    springs=0,
    dampers=0, 
    limits=false,
    joint_limits=Dict(),
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact=true,
    contact_radius=0,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}()
    block = Box(edge_length, edge_length, edge_length, mass; color, name=:block)
    bodies = [block]

    joint = JointConstraint(PlanarAxis(origin, block, X_AXIS); name=:joint)
    joints = [joint]

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
    origin = mechanism.origin
    bodies = mechanism.bodies
    joints = mechanism.joints
    contacts = ContactConstraint{T}[]

    if contact
        # Corner vectors
        corners = [
            [[0,  edge_length / 2,  edge_length / 2]]
            [[0,  edge_length / 2, -edge_length / 2]]
            [[0, -edge_length / 2,  edge_length / 2]]
            [[0, -edge_length / 2, -edge_length / 2]]
        ]

        n = length(corners)
        normal = [Z_AXIS for i = 1:n]
        contact_radius = [contact_radius for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts = contact_constraint(block, normal;
            friction_coefficient,
            contact_origins=corners,
            contact_radius,
            contact_type,
            names=[Symbol(:contact, i) for i=1:8])
    end

    mechanism = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_block2d!(mechanism)
        
    # construction finished
    return mechanism
end

function initialize_block2d!(mechanism::Mechanism;
    position = zeros(2), orientation=0, velocity=zeros(2), angular_velocity=0)

    zero_velocity!(mechanism)
    zero_coordinates!(mechanism)
    
    body = mechanism.bodies[1]
    edge_length = body.shape.xyz[3]/2

    if length(mechanism.contacts) > 0
        model = mechanism.contacts[1].model
        offset = model.collision.contact_radius
    else
        offset = 0.0
    end

    height = edge_length + offset

    set_maximal_configurations!(body,
        x=[0;position] + Z_AXIS*height, q=RotX(orientation))
    set_maximal_velocities!(body,
        v=[0;velocity], ω=[angular_velocity;0;0])

    return
end