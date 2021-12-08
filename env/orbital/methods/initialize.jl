function getorbital(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0, Nlink::Int = 5) where {T}
    # Parameters
    ex = [0; 0; 1]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()

    links = [Box(r, r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Fixed(origin, links[1]; p2 = vert11))
    if Nlink > 1
        eqcs = [
            jointb1;
            [EqualityConstraint(Orbital(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink]
            ]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function initializeorbital!(mechanism::Mechanism; ϕx::T = pi/4, ϕy::T = pi/8) where {T}

    link1 = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = vert11, Δq = UnitQuaternion(RotX(0.0)))

    previd = link1.id
    setPosition!(mechanism, collect(mechanism.eqconstraints)[2], [ϕx, ϕy])
end