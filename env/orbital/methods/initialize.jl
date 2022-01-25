function getorbital(; Δt::T = 0.01, g::T = -9.81, spring = 0.0, damper = 0.0, Nb::Int = 5) where {T}
    # Parameters
    ex = [0; 0; 1]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()

    bodies = [Box(r, r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Fixed(origin, bodies[1]; p2 = vert11))
    if Nb > 1
        eqcs = [
            jointb1;
            [JointConstraint(Orbital(bodies[i - 1], bodies[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nb]
            ]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end

function initializeorbital!(mechanism::Mechanism; ϕx::T = pi/4, ϕy::T = pi/8) where {T}
    body1 = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, body1, p2 = vert11, Δq = UnitQuaternion(RotX(0.0)))

    previd = body1.id
    setPosition!(mechanism, collect(mechanism.eqconstraints)[2], [ϕx, ϕy])

    zeroVelocity!(mechanism)
end
