function getsphere(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, radius=0.5,
        contact::Bool=true, contact_type::Symbol=:contact) where T
    origin = Origin{T}(name=:origin)
    mass = 1.0
    bodies = [Sphere(radius, mass, name=:sphere)]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = :floating_joint)]
    mechanism = Mechanism(origin, bodies, eqcs, Δt = Δt, g = g)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        ineqcs = [contactconstraint(getbody(mechanism, :sphere), normal, cf=cf,
            p=contact, offset=[0,0,radius], contact_type=contact_type)]
        setPosition!(mechanism, geteqconstraint(mechanism, :floating_joint), [0;0;radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, eqcs, ineqcs, g=g, Δt=Δt)
    end
    return mechanism
end

function initializesphere!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}
    r = collect(mechanism.bodies)[1].shape.r
    eqc = geteqconstraint(mechanism, :floating_joint)
    zeroVelocity!(mechanism)
    setPosition!(mechanism, eqc, [x+[0,0,r] rotation_vector(q)])
    setVelocity!(mechanism, eqc, [v; ω])
end
