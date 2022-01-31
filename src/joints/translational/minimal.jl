@inline function minimal_coordinates(joint::Translational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * position_error(joint, xa, qa, xb, qb)
end

@inline function minimal_velocities(joint::Translational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.x2[1], statea.q2[1], statea.v15, statea.ϕ15, stateb.x2[1], stateb.q2[1], stateb.v15, stateb.ϕ15)
end

@inline function minimal_velocities(joint::Translational, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
    Δv = vrotate(Δvw, inv(qa)) # in the a frame
    return nullspace_mask(joint) * Δv
end
