@inline function getPositionDelta(joint::Rotational, body1::Node, body2::Node, θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspacemat(joint)) * θ
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end

    Δq = q * joint.qoffset # in body1 frame
    return Δq
end

@inline function getVelocityDelta(joint::Rotational, body1::Node, body2::Node, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = ω # in body1 frame
    return Δω
end

@inline function minimalCoordinates(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimalCoordinates(joint, statea.q2[1], stateb.q2[1])
end

@inline function minimalCoordinates(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    q = qa \ qb / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end

@inline function minimalVelocities(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimalVelocities(joint, statea.q2[1], statea.ϕ15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimalVelocities(joint::Rotational, qa::UnitQuaternion,
        ϕa::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    return nullspacemat(joint) * (vrotate(ϕb, qa \ qb) - ϕa) # in body1's frame
end

