function gettippetop(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, contact::Bool=true, contact_type::Symbol=:contact) where T
    origin = Origin{T}(name="origin")
    radius = 0.5
    mass = 1.0
    α = 0.2
    bodies = [Sphere(radius, mass, name="sphere1"), Sphere(radius*α, mass*α^3, name="sphere2")]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = "floating_joint"),
        EqualityConstraint(Fixed(bodies[1], bodies[2], p1=[0,0,radius], p2=zeros(3)), name = "fixed_joint"),]
    mechanism = Mechanism(origin, bodies, eqcs, Δt = Δt, g = g)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        ineqcs = [
            contactconstraint(getbody(mechanism, "sphere1"), normal, cf=cf, p=contact, offset=[0,0,radius], contact_type=contact_type),
            contactconstraint(getbody(mechanism, "sphere2"), normal, cf=cf, p=contact, offset=[0,0,radius*α], contact_type=contact_type)
            ]
        setPosition!(mechanism, geteqconstraint(mechanism, "floating_joint"), [0;0;radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, eqcs, ineqcs, g=g, Δt=Δt)
    end
    return mechanism
end

function initializetippetop!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}

    eqc2 = geteqconstraint(mechanism, "fixed_joint")
    radius = eqc2.constraints[1].vertices[1][3]
    origin = mechanism.origin
    body1 = getbody(mech, "sphere1")
    body2 = getbody(mech, "sphere2")

    zeroVelocity!(mechanism)
    # setPosition!(mechanism, eqc, [x; rotation_vector(q)])
    # setVelocity!(mechanism, eqc, [v; ω])
    setPosition!(origin, body1; p1 = [0;0;radius], p2 = [0;0;0], Δx = x, Δq = q)
    setPosition!(body1,  body2; p1 = [0;0;radius], p2 = [0;0;0], Δx = [0;0;0], Δq = one(UnitQuaternion))
    setVelocity!(origin, body1; p1 = [0;0;radius], p2 = [0;0;0], Δv = v, Δω = ω)
    setVelocity!(body1,  body2; p1 = [0;0;radius], p2 = [0;0;0], Δv = [0;0;0], Δω = [0;0;0])
    return nothing
end
