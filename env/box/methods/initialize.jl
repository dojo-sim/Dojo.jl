function getbox(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, radius = 0.0, side = 0.5,
    contact::Bool = true,
    conetype = :soc,
    # conetype = :impact,
    # conetype = :linear,
    color=RGBA(1.0, 0.0, 0.0, 1.0),
    mode = :box)  where {T}
    # Parameters
    origin = Origin{T}()
    link1 = Box(side, side, side, 1., color=color)
    joint0to1 = EqualityConstraint(Floating(origin, link1))
    links = [link1]
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

        if conetype == :soc
            contineqcs = contactconstraint(link1, normal, cf, p = corners, offset=offset)
            mech = Mechanism(origin, links, eqcs, contineqcs, g = g, Δt = Δt)
        elseif conetype == :impact
            impineqcs = impactconstraint(link1, normal, p = corners, offset=offset)
            mech = Mechanism(origin, links, eqcs, impineqcs, g = g, Δt = Δt)
        elseif conetype == :linear
            linineqcs = linearcontactconstraint(link1, normal, cf, p=corners, offset=offset)
            mech = Mechanism(origin, links, eqcs, linineqcs, g = g, Δt = Δt)
        else
            error("Unknown conetype")
        end
    else
        mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    end
    return mech
end

function initializebox!(mechanism::Mechanism;
        x::AbstractVector{T} = [0,0,1.],
        q::UnitQuaternion{T} = UnitQuaternion(1.,0,0,0),
        v::AbstractVector{T} = [1,.3,.2],
        ω::AbstractVector{T} = [2.5,-1,2]) where {T}
    halfside = mechanism.bodies.values[1].shape.xyz[1]/2
    if length(mechanism.ineqconstraints.values) > 0
        bound = mechanism.ineqconstraints.values[1].constraints[1]
        offset = bound.offset[3]
    else
        offset = 0.0
    end
    z = halfside + offset
    body = mechanism.bodies.values[1]
    setPosition!(body, x = x + [0,0,z], q = q)
    setVelocity!(body, v = v, ω = ω)
end
