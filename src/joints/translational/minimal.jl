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
    return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15, stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimal_velocities(joint::Translational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
    Δv = vrotate(Δvw, inv(qa)) # in the a frame
    return nullspace_mask(joint) * Δv
end

function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol, joint::Translational,
	xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
	jacobian_relative==:parent && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_coordinates(joint, x, qa, xb, qb), xa) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, xa, UnitQuaternion(q..., false), xb, qb), vector(qa))]) # * LVᵀmat(qa)]) 
	jacobian_relative==:child && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_coordinates(joint, xa, qa, x, qb), xb) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, xa, qa, xb, UnitQuaternion(q..., false)), vector(qb))]) # * LVᵀmat(qb)]) 
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
	joint::Translational, xa::AbstractVector, va::AbstractVector,
	qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
	vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)

	(jacobian_relative == :parent) && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_velocities(joint, x, va, qa, ωa, xb, vb, qb, ωb), xa) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, xa, va, UnitQuaternion(q..., false), ωa, xb, vb, qb, ωb), vector(qa))]) # * LVᵀmat(qa)])
	(jacobian_relative == :child) && (return [FiniteDiff.finite_difference_jacobian(x -> minimal_velocities(joint, xa, va, qb, ωa, x, vb, qb, ωb), xb) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, xa, va, qa, ωa, xb, vb, UnitQuaternion(q..., false), ωb), vector(qb))]) # * LVᵀmat(qb)])
	return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
	joint::Translational, xa::AbstractVector, va::AbstractVector,
	qa::UnitQuaternion, ωa::AbstractVector, xb::AbstractVector,
	vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	(jacobian_relative == :parent) && (return [FiniteDiff.finite_difference_jacobian(v -> minimal_velocities(joint, xa, v, qa, ωa, xb, vb, qb, ωb), va) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, xa, va, qa, ω, xb, vb, qb, ωb), ωa)])
	(jacobian_relative == :child) && (return [FiniteDiff.finite_difference_jacobian(v -> minimal_velocities(joint, xa, va, qa, ωa, xb, v, qb, ωb), vb) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ω), ωb)])
	return
end