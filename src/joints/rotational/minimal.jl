@inline function get_position_delta(joint::Rotational, body1::Node, body2::Node, θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspace_mask(joint)) * θ
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end

    Δq = q * joint.qoffset # in body1 frame
    return Δq
end

@inline function get_velocity_delta(joint::Rotational, body1::Node, body2::Node, ω::SVector)
    ω = zerodimstaticadjoint(nullspace_mask(joint)) * ω
    Δω = ω # in body1 frame
    return Δω
end

@inline function minimal_coordinates(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.q2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    q = qa \ qb / joint.qoffset
    return nullspace_mask(joint) * rotation_vector(q)
end

@inline function minimal_velocities(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.q2[1], statea.ϕ15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimal_velocities(joint::Rotational, qa::UnitQuaternion,
        ϕa::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    return nullspace_mask(joint) * (vrotate(ϕb, qa \ qb) - ϕa) # in body1's frame
end

