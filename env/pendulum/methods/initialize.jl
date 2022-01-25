function getpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
        spring = 0.0, damper = 0.0, spring_offset = szeros(1)) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    body1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_body1 = EqualityConstraint(Revolute(origin, body1,
        joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset))
    bodies = [body1]
    eqcs = [joint_between_origin_and_body1]

    mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end

function initializependulum!(mechanism::Mechanism; ϕ1::T = 0.7, ω1::T = 0.0) where {T}
    body = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    p2 = eqc.constraints[1].vertices[2]
    p1 = eqc.constraints[1].vertices[1]
    q1 = UnitQuaternion(RotX(ϕ1))
    setPosition!(mechanism.origin, body, p1 = p1, p2 = p2, Δq = q1)
    setVelocity!(mechanism.origin, body, p1 = p1, p2 = p2, Δω = [ω1,0,0])
end

function getnpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
        spring::T = 0.0, damper::T = 0.0, Nb::Int = 5,
        basetype::Symbol = :Revolute, jointtype::Symbol = :Revolute) where {T}
    # Parameters
    ex = [1.; 0; 0]
    r = 0.05
    vert11 = [0; 0; l/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    bodies = [Box(r, r, l, m, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = EqualityConstraint(Prototype(basetype, origin, bodies[1], ex; p2 = vert11, spring = spring, damper = damper))
    if Nb > 1
        eqcs = [EqualityConstraint(Prototype(jointtype, bodies[i - 1], bodies[i], ex; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nb]
        eqcs = [jointb1; eqcs]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt)
    return mech
end

function initializenpendulum!(mechanism::Mechanism; ϕ1::T = pi/4, ω = [0.0, 0.0, 0.0],
    Δv::AbstractVector{T} = [0, 0, 0.], Δω::AbstractVector{T} = [0, 0, 0.]) where {T}

    body1 = mechanism.bodies[1]
    eqc = mechanism.eqconstraints[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, body1, p2 = vert11, Δq = UnitQuaternion(RotX(ϕ1)))
    setVelocity!(body1, ω = ω)

    previd = body1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = 1/i*Δω)
        previd = body.id
    end
end
