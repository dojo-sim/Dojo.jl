function get_block(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.0,
    side=0.5,
    contact=true,
    contact_type=:nonlinear,
    color=RGBA(0.3, 0.3, 0.3, 1.0),
    mode=:box,
    T=Float64)

    # Parameters
    origin = Origin{T}()
    pbody = Box(side, side, side, 1.0,
        color=color)
    joint0to1 = JointConstraint(Floating(origin, pbody))
    bodies = [pbody]
    joints = [joint0to1]

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

        contacts = contact_constraint(pbody, normal;
            friction_coefficient,
            contact_origins=corners,
            contact_radius,
            contact_type)

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
        x=[0.0, 0.0, 1.0],
        q=Quaternion(1.0, 0.0, 0.0, 0.0),
        v=[1.0, 0.3, 0.2],
        ω=[2.5, -1.0, 2.0]) where T

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
        x=x + [0.0, 0.0, z],
        q=q)
    set_maximal_velocities!(body,
        v=v,
        ω=ω)
end
