function get_particle(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=0.2,
    side=0.5,
    contact=true,
    contact_type=:nonlinear,
    T=Float64)

    # Parameters
    axis = [1.0, 0.0, 0.0]

    origin = Origin{T}()
    pbody = Box(side, side, side, 1.,
        color=RGBA(1., 1., 0.), name=:particle)
    joint = JointConstraint(FixedOrientation(origin, pbody), name=:fixed_orientation)
    bodies = [pbody]
    joints = [joint]

    if contact
        # Corner vectors
        corners = [[0.0, 0.0, -side / 2.0]]
        n = length(corners)
        normal = [[0.0, 0.0, 1.0] for i = 1:n]
        contact_radius = [0.0 for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts = contact_constraint(pbody, normal,
            friction_coefficient=friction_coefficient,
            contact_origins=corners,
            contact_radius=contact_radius,
            contact_type=contact_type)

        mech = Mechanism(origin, bodies, joints, contacts,
            gravity=gravity,
            timestep=timestep)
    else
        mech = Mechanism(origin, bodies, joints,
            gravity=gravity,
            timestep=timestep)
    end
    return mech
end

function initialize_particle!(mechanism::Mechanism{T};
    position=[0.0, 0.0, 1.0],
    velocity=[0.0, 0.0, 0.0]) where T

    position += [0,0,0.25]
    body = mechanism.bodies[1]
    set_maximal_configurations!(body, x=position)
    set_maximal_velocities!(body, v=velocity)
end
