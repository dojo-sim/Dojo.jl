function getsphere(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, contact::Bool=true) where T
    origin = Origin{T}(name="origin")
    radius = 0.5
    mass = 1.0
    bodies = [Sphere(radius, mass, name="sphere")]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = "floating_joint")]
    mechanism = Mechanism(orig, bodies, eqcs, Δt = Δt, g = g)

    if contact
        contact = [0,0,0.0]
        normal = [0,0,1.0]
        contineqcs = [contactconstraint(getbody(mechanism, "sphere"), normal, cf, p=contact)]
        setPosition!(mechanism, geteqconstraint(mechanism, "floating_joint"), [0;0;radius;zeros(3)])
        mechanism = Mechanism(origin, bodies, eqcs, contineqcs, g=g, Δt=Δt)
    end
    return mechanism
end

function initializesphere!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}

    eqc = geteqconstraint(mechanism, "floating_joint")
    zeroVelocity!(mechanism)
    setPosition!(mechanism, eqc, [x; rotation_vector(q)])
    setVelocity!(mechanism, eqc, [v; ω])
end
