function gettennisracket(; Δt::T=0.01, g::T=-9.81) where T
    origin = Origin{T}(name="origin")
    mass = 1.0
    r = 0.1
    h = 1.0
    bodies = [Box(h/25, h/2, h, mass, color = RGBA(1., 0., 0.), name = :box)]
    eqcs = [EqualityConstraint(Floating(origin, bodies[1]), name = :floating_joint)]
    mechanism = Mechanism(origin, bodies, eqcs, Δt = Δt, g = g)
    return mechanism
end

function initializetennisracket!(mechanism::Mechanism; x::AbstractVector{T}=zeros(3),
        q::UnitQuaternion{T}=one(UnitQuaternion), v::AbstractVector{T}=zeros(3),
        ω::AbstractVector{T}=zeros(3)) where {T}

    eqc = geteqconstraint(mechanism, :floating_joint)
    zeroVelocity!(mechanism)
    setPosition!(mechanism, eqc, [x; rotation_vector(q)])
    setVelocity!(mechanism, eqc, [v; ω])
end
