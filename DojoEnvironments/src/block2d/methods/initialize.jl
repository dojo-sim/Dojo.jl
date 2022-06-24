function get_block2d(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.8,
    radius=0.0,
    side=0.5,
    contact=true,
    contact_type=:nonlinear,
    mode=:block2d,
    T=Float64)

    # Parameters
    axis = [1.0, 0.0, 0.0]

    origin = Origin{T}()
    block = Box(side, side, side, 1.; color=RGBA(1., 1., 0.), name=:block)
    joint1 = JointConstraint(PlanarAxis(origin, block, axis), name=:joint)
    bodies = [block]
    joints = [joint1]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0.0, 0.0, 0.0]]
        elseif mode == :block2d
            corners = [
                [[0.0,  side / 2.0,  side / 2.0]]
                [[0.0,  side / 2.0, -side / 2.0]]
                [[0.0, -side / 2.0,  side / 2.0]]
                [[0.0, -side / 2.0, -side / 2.0]]
            ]
        else
            @error "incorrect mode specified, try :particle or :block2d"
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

function initialize_block2d!(mechanism::Mechanism{T};
    position=[0.0, 1.0],
    velocity=[0.0, 0.0],
    orientation=0.0,
    angular_velocity=0.0) where T

    if length(mechanism.contacts) > 0
        model = mechanism.contacts[1].model
        side = model.collision.contact_origin[2]
        offset = model.collision.contact_radius
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
