function getsnake(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, contact::Bool=true,
    contact_type=:contact, spring=0.0, damper=0.0, Nb::Int=2,
    jointtype::Symbol=:Spherical, h::T=1.0, r::T=0.05) where {T}

    # Parameters
    ex = [1.;0.;0.]
    ey = [0.;1.;0.]
    ez = [0.;0.;1.]

    vert11 = [0.;0.;h / 2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    # bodies = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]
    bodies = [Box(3r, 2r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Floating(origin, bodies[1], spring = 0.0, damper = 0.0))
    if Nb > 1
        eqcs = [JointConstraint(Prototype(jointtype, bodies[i - 1], bodies[i], ex; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nb]
        eqcs = [jointb1; eqcs]
    else
        eqcs = [jointb1]
    end

    if contact
        n = Nb
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(n)

        ineqcs1 = contactconstraint(bodies, normal, cf=cf, p = fill(vert11, n), contact_type=contact_type) # we need to duplicate point for prismatic joint for instance
        ineqcs2 = contactconstraint(bodies, normal, cf=cf, p = fill(vert12, n), contact_type=contact_type)
        mech = Mechanism(origin, bodies, eqcs, [ineqcs1; ineqcs2], g=g, Δt=Δt, spring=spring, damper=damper)
    else
        mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializesnake!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T}=[0,-0.5,0],
    v::AbstractVector{T}=zeros(3), ω::AbstractVector{T}=zeros(3),
    Δω::AbstractVector{T}=zeros(3), Δv::AbstractVector{T}=zeros(3),
    q1::UnitQuaternion{T}=UnitQuaternion(RotX(0.6 * π))) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    body1 = bodies[1]
    # h = body1.shape.rh[2]
    h = body1.shape.xyz[3]
    vert11 = [0.;0.; h/2]
    vert12 = -vert11
    # set position and velocities
    setPosition!(mechanism.origin, body1, p2 = x, Δq = q1)
    setVelocity!(body1, v = v, ω = ω)

    previd = body1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = Δω)
        previd = body.id
    end
end
