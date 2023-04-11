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
    contacts = ContactConstraint{T}[]

    if contact
        names = [Symbol(:contact, i) for i=1:4]
        n = length(names)
        normals = fill(Z_AXIS,n)
        friction_coefficients = fill(friction_coefficient,n)
        contact_origins = [
            [0,  edge_length / 2,  edge_length / 2],
            [0,  edge_length / 2, -edge_length / 2],
            [0, -edge_length / 2,  edge_length / 2],
            [0, -edge_length / 2, -edge_length / 2],
        ]
        contact_radii = fill(contact_radius,n)

        contacts = [contacts;contact_constraint(block, normals; friction_coefficients, contact_origins, contact_radii, contact_type, names)]
    end

    mechanism = Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, contacts;
        gravity, timestep, input_scaling)

    # zero configuration
    initialize_block2d!(mechanism)
        
    # construction finished
    return mechanism
end

function initialize_block2d!(mechanism::Mechanism;
    position = [0;1], orientation=0, velocity=zeros(2), angular_velocity=randn())

    zero_velocities!(mechanism)
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
        x=[0;position] + Z_AXIS*height, q=Dojo.RotX(orientation))
    set_maximal_velocities!(body,
        v=[0;velocity], ω=[angular_velocity;0;0])

    return
end