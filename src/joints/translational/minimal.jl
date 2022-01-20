## Minimal coordinate calculation
@inline function minimalCoordinates(joint::Translational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalCoordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimalCoordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspacemat(joint) * position_error(joint, xa, qa, xb, qb)
end

@inline function minimalVelocities(joint::Translational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalVelocities(joint, statea.x2[1], statea.q2[1], statea.v15, statea.ϕ15, stateb.x2[1], stateb.q2[1], stateb.v15, stateb.ϕ15)
end

@inline function minimalVelocities(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa))
    Δv = vrotate(Δvw, inv(qa))
    return nullspacemat(joint) * Δv
end

