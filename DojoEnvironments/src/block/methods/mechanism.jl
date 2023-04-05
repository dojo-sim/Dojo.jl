function get_block(;
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    mass=1,
    edge_length=0.5,
    color=RGBA(0.9, 0.9, 0.9, 1),
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

    joint = JointConstraint(Floating(origin, block))
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
            [[ edge_length / 2;  edge_length / 2; -edge_length / 2]]
            [[ edge_length / 2; -edge_length / 2; -edge_length / 2]]
            [[-edge_length / 2;  edge_length / 2; -edge_length / 2]]
            [[-edge_length / 2; -edge_length / 2; -edge_length / 2]]
            [[ edge_length / 2;  edge_length / 2;  edge_length / 2]]
            [[ edge_length / 2; -edge_length / 2;  edge_length / 2]]
            [[-edge_length / 2;  edge_length / 2;  edge_length / 2]]
            [[-edge_length / 2; -edge_length / 2;  edge_length / 2]]
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
    zero_coordinates!(mechanism)
    set_minimal_coordinates!(mechanism, joint, [0; 0; edge_length/2+contact_radius[1]; 0; 0; 0])
    
    # construction finished
    return mechanism
end

function initialize_block!(mechanism::Mechanism{T};
    position=Z_AXIS,
    orientation=Quaternion(1, 0, 0, 0),
    velocity=[1, 0.3, 0.2],
    angular_velocity=[2.5, -1, 2]) where T

body = mechanism.bodies[1]

halfside = body.shape.xyz[1] / 2.0

if length(mechanism.contacts) > 0
    model = mechanism.contacts[1].model
    offset = model.collision.contact_radius
else
    offset = 0.0
end

z = halfside + offset

set_maximal_configurations!(body,
    x=position + [0, 0, z],
    q=orientation)
set_maximal_velocities!(body,
    v=velocity,
    ω=angular_velocity)
end
