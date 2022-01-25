function getslider(; Δt::T = 0.01, g::T = -9.81, spring = 0.0, damper = 0.0) where {T}
    # Parameters
    joint_axis = [0; 0; 1.0]
    length1 = 1.0
    width, depth = 0.1, 0.1
    p2 = [0; 0; length1/2] # joint connection point

    # Links
    origin = Origin{T}()
    body1 = Box(width, depth, length1, length1)

    # Constraints
    joint_between_origin_and_body1 = JointConstraint(Prismatic(origin, body1, joint_axis; p2=p2, spring = spring, damper = damper))
    bodies = [body1]
    eqcs = [joint_between_origin_and_body1]

    mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end

function initializeslider!(mechanism::Mechanism; z1::T = 0.0) where {T}
    body = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.joints)[1]
    p2 = eqc.constraints[1].vertices[2]
    set_position(mechanism.origin, body, p2 = p2 - [0, 0, z1])
end

function getnslider(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0, Nb::Int = 5) where {T}
    # Parameters
    ex = [0; 0; 1.0]
    h = 1.
    r = .05
    vert11 = [0; r; 0.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Prismatic(origin, bodies[1], ex; p2 = 0*vert11))
    if Nb > 1
        eqcs = [
            jointb1;
            [JointConstraint(Prismatic(bodies[i - 1], bodies[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nb]
            ]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt)
    return mech
end

function initializenslider!(mechanism::Mechanism; z1::T = 0.2, Δz = 0.0) where {T}
    body1 = collect(mechanism.bodies)[1]
    # set position and velocities
    set_position(mechanism.origin, body1, p1 = [0, 0, z1])

    previd = body1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        set_position(get_body(mechanism, previd), body, p1 = [0, -0.1, Δz])
        previd = body.id
    end

    zeroVelocity!(mechanism)
end
