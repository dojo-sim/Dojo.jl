function getdzhanibekov(; Δt::T = 0.01, g::T = -9.81, h::T = 1.0, r::T = 0.05) where {T}
    # Parameters
    p1 = [3r/2,0,0.]
    p2 = [0,0,-h/6]
    qoffset = UnitQuaternion(RotY(π/2))

    # Links
    origin = Origin{T}()
    link1 = Box(3r, 2r, h, h, color = RGBA(1., 0., 0.))
    link2 = Box(3r/3, 2r/3, h/3, h/3, color = RGBA(1., 0., 0.))
    links = [link1, link2]

    # Constraints
    eqc1 = EqualityConstraint(Floating(origin, links[1]))
    eqc2 = EqualityConstraint(Fixed(link1, link2; p1 = p1, p2 = p2, qoffset = qoffset))
    eqcs = [eqc1, eqc2]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function initializedzhanibekov!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T} = [0,0,0.],
    v::AbstractVector{T} = zeros(3), ω::AbstractVector{T} = zeros(3),
    q::UnitQuaternion{T} = UnitQuaternion(RotX(0.0 * π))) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    link1 = bodies[1]
    link2 = bodies[2]
    # h = link1.shape.rh[2]
    h = link1.shape.xyz[3]
    r = link1.shape.xyz[2]/2
    p1 = [3r/2,0,0.]
    p2 = [0,0,-h/6]
    qoffset = UnitQuaternion(RotY(π/2))

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = x, Δq = q)
    setVelocity!(link1, v = v, ω = ω)
    setPosition!(link1, link2, p1 = p1, p2 = p2, Δq = qoffset)
    setVelocity!(link1, link2, p1 = p1, p2 = p2)
end