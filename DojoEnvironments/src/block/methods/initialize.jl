function get_block(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.0,
    side=0.5,
    contact=true,
    contact_type=:nonlinear,
    color=RGBA(0.9, 0.9, 0.9, 1.0),
    mode=:box,
    T=Float64)

    # Parameters
    origin = Origin{T}()
    block = Box(side, side, side, 1.0, color=color, name=:block)
    bodies = [block]
    joints = [JointConstraint(Floating(origin, block))]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0.0, 0.0, 0.0]]
        elseif mode == :box
            corners = [
                [[ side / 2.0;  side / 2.0; -side / 2.0]]
                [[ side / 2.0; -side / 2.0; -side / 2.0]]
                [[-side / 2.0;  side / 2.0; -side / 2.0]]
                [[-side / 2.0; -side / 2.0; -side / 2.0]]
                [[ side / 2.0;  side / 2.0;  side / 2.0]]
                [[ side / 2.0; -side / 2.0;  side / 2.0]]
                [[-side / 2.0;  side / 2.0;  side / 2.0]]
                [[-side / 2.0; -side / 2.0;  side / 2.0]]
            ]
        else
            @error "incorrect mode specified, try :particle or :box"
        end
        n = length(corners)
        normal = [[0.0, 0.0, 1.0] for i = 1:n]
        contact_radius = [radius for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts = contact_constraint(block, normal;
            friction_coefficient,
            contact_origins=corners,
            contact_radius,
            contact_type,
            names=[Symbol(:contact, i) for i=1:8])

        mech = Mechanism(origin, bodies, joints, contacts;
            gravity,
            timestep)
    else
        mech = Mechanism(origin, bodies, joints;
            gravity,
            timestep)
    end
    return mech
end

function initialize_block!(mechanism::Mechanism{T};
        position=[0.0, 0.0, 1.0],
        orientation=Quaternion(1.0, 0.0, 0.0, 0.0),
        velocity=[1.0, 0.3, 0.2],
        angular_velocity=[2.5, -1.0, 2.0]) where T

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
        x=position + [0.0, 0.0, z],
        q=orientation)
    set_maximal_velocities!(body,
        v=velocity,
        ω=angular_velocity)
end
