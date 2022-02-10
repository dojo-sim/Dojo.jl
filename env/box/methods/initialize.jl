function get_box(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=0.8, radius=0.0, side=0.5,
    contact::Bool=true,
    contact_type=:contact,
    # contact_type=:linear_contact,
    # contact_type=:impact,
    color=RGBA(0.0, 0.0, 0.0, 1.0),
    mode=:box)  where T
    # Parameters
    origin = Origin{T}()
    body1 = Box(side, side, side, 1., color=color)
    joint0to1 = JointConstraint(Floating(origin, body1))
    bodies = [body1]
    joints = [joint0to1]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0,0,0.]]
        elseif mode == :box
            corners = [
                [[ side/2;  side/2; -side/2]]
                [[ side/2; -side/2; -side/2]]
                [[-side/2;  side/2; -side/2]]
                [[-side/2; -side/2; -side/2]]
                [[ side/2;  side/2;  side/2]]
                [[ side/2; -side/2;  side/2]]
                [[-side/2;  side/2;  side/2]]
                [[-side/2; -side/2;  side/2]]
            ]
        else
            @error "incorrect mode specified, try :particle or :box"
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

function initialize_box!(mechanism::Mechanism;
        x::AbstractVector{T} = [0,0,1.],
        q::UnitQuaternion{T} = UnitQuaternion(1.,0,0,0),
        v::AbstractVector{T} = [1,.3,.2],
        ω::AbstractVector{T} = [2.5,-1,2]) where T
        
    body = mechanism.bodies[1]

    halfside = body.shape.xyz[1] / 2

    if length(mechanism.contacts) > 0
        model = mechanism.contacts[1].model
        offset = model.offset[3]
    else
        offset = 0.0
    end

    z = halfside + offset
    set_position!(body, x = x + [0,0,z], q = q)
    set_velocity!(body, v = v, ω = ω)
end
