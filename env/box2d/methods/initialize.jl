function get_box2d(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=0.8, radius=0.0, side=0.5,
    contact::Bool=true,
    contact_type=:contact,
    mode=:box2d)  where T
    # Parameters
    axis = [1,0,0.]

    origin = Origin{T}()
    body1 = Box(side, side, side, 1., color = RGBA(1., 1., 0.))
    joint1 = JointConstraint(PlanarAxis(origin, body1, axis))
    bodies = [body1]
    joints = [joint1]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0,0,0.]]
        elseif mode == :box2d
            corners = [
                [[0,  side/2,  side/2]]
                [[0,  side/2, -side/2]]
                [[0, -side/2,  side/2]]
                [[0, -side/2, -side/2]]
            ]
        else
            @error "incorrect mode specified, try :particle or :box2d"
        end
        n = length(corners)
        normal = [[0,0,1.0] for i = 1:n]
        offset = [[0,0,radius] for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts = contact_constraint(body1, normal, friction_coefficient=friction_coefficient, contact_points=corners, offset=offset, contact_type=contact_type)
        mech = Mechanism(origin, bodies, joints, contacts, gravity=gravity, timestep=timestep)
    else
        mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep)
    end
    return mech
end

function initialize_box2d!(mechanism::Mechanism{T}; x::AbstractVector=[0,1.],
    v::AbstractVector=[0,0], θ::T=0.0, ω::T=0.0) where T
    if length(mechanism.contacts) > 0
        bound = mechanism.contacts[1].constraints[1]
        side = bound.p[2]
        offset = bound.offset[3]
        z = side + offset
    else
        z = 0.0
    end
    body = mechanism.bodies[1]
    set_position!(body, x = [0;x] + [0,0,z], q = UnitQuaternion(RotX(θ)))
    set_velocity!(body, v = [0;v], ω = [ω,0,0])
end
