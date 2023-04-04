function get_block2d(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1
    edge_length=0.5,
    color=RGBA(1, 1, 0.),
    springs=0.0,
    dampers=0.0, 
    limits=false,
    joint_limits=Dict(),
    friction_coefficient=0.8,
    contact=true,
    contact_radius=0.0,
    contact_type=:nonlinear,
    T=Float64)

    # mechanism
    origin = Origin{T}()
    block = Box(edge_length, edge_length, edge_length, mass; color, name=:block)
    bodies = [block]

    joint = JointConstraint(PlanarAxis(origin, block, X_AXIS), name=:joint)
    joints = [joint]

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
        # Corner vectors
        corners = [
            [[0.0,  edge_length / 2.0,  edge_length / 2.0]]
            [[0.0,  edge_length / 2.0, -edge_length / 2.0]]
            [[0.0, -edge_length / 2.0,  edge_length / 2.0]]
            [[0.0, -edge_length / 2.0, -edge_length / 2.0]]
        ]

        n = length(corners)
        normal = [Z_AXIS for i = 1:n]
        contact_radii = [contact_radius for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts = contact_constraint(block, normal;
            friction_coefficient,
            contact_origins=corners,
            contact_radii,
            contact_type,
            names=[Symbol(:contact, i) for i=1:8])
    end

    mech = Mechanism(origin, bodies, joints, contacts;
        gravity, timestep, input_scaling)
        
    # construction finished
    return mech
end

function initialize_block2d!(mechanism::Mechanism{T};
    position=[0.0, 1.0],
    velocity=[0.0, 0.0],
    orientation=0.0,
    angular_velocity=0.0) where T

    if length(mechanism.contacts) > 0
        model = mechanism.contacts[1].model
        side = model.collision.contact_origin[2]
        offset = model.collision.contact_radii
        z = side + offset
    else
        z = 0.0
    end

    body = mechanism.bodies[1]

    set_maximal_configurations!(body,
        x=[0.0; position] + [0.0, 0.0 , z],
        q=RotX(orientation))
    set_maximal_velocities!(body,
        v=[0.0; velocity],
        ω=[angular_velocity, 0.0, 0.0])
end
