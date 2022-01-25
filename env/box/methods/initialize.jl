function getbox(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, radius=0.0, side=0.5,
    contact::Bool=true,
    contact_type=:contact,
    # contact_type=:linear_contact,
    # contact_type=:impact,
    color=RGBA(1.0, 0.0, 0.0, 1.0),
    mode=:box)  where {T}
    # Parameters
    origin = Origin{T}()
    body1 = Box(side, side, side, 1., color=color)
    joint0to1 = JointConstraint(Floating(origin, body1))
    bodies = [body1]
    eqcs = [joint0to1]

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
        cf = cf * ones(n)

        ineqcs = contactconstraint(body1, normal, cf=cf, p=corners, offset=offset, contact_type=contact_type)
        mech = Mechanism(origin, bodies, eqcs, ineqcs, g=g, Δt=Δt)
    else
        mech = Mechanism(origin, bodies, eqcs, g=g, Δt=Δt)
    end
    return mech
end

function initializebox!(mechanism::Mechanism;
        x::AbstractVector{T} = [0,0,1.],
        q::UnitQuaternion{T} = UnitQuaternion(1.,0,0,0),
        v::AbstractVector{T} = [1,.3,.2],
        ω::AbstractVector{T} = [2.5,-1,2]) where {T}
        
    body = mechanism.bodies[1]

    halfside = body.shape.xyz[1] / 2

    if length(mechanism.ineqconstraints) > 0
        bound = mechanism.ineqconstraints[1].constraints[1]
        offset = bound.offset[3]
    else
        offset = 0.0
    end

    z = halfside + offset
    setPosition!(body, x = x + [0,0,z], q = q)
    setVelocity!(body, v = v, ω = ω)
end
