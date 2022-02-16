################################################################################
# Displacements
################################################################################
@inline function displacement(joint::Translational, xa::AbstractVector,
		qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    rotate && (return vrotate(d, inv(qa))) : (return d)
end

@inline function displacement_jacobian_configuration(relative::Symbol, joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac=true) where T

    vertices = joint.vertices

    if relative == :parent
        d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
        X = -rotation_matrix(inv(qa))
        Q = -rotation_matrix(inv(qa)) * ∂qrotation_matrix(qa, vertices[1])
        Q += ∂qrotation_matrix_inv(qa, d)
        attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
        X = rotation_matrix(inv(qa))
        Q = rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, vertices[2])
        attjac && (Q *= LVᵀmat(qb))
    end

    return X, Q
end

################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * displacement(joint, xa, qa, xb, qb)
end

@inline function minimal_coordinates_jacobian_configuration(relative::Symbol, joint::Translational,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac::Bool=true)
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=attjac)
	return nullspace_mask(joint) * [X Q]
end

################################################################################
# Velocities
################################################################################
@inline function minimal_velocities(joint::Translational,
		xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
		xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
		timestep)
	A = nullspace_mask(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	# Coordinates
	Δx = A * displacement(joint, xa, qa, xb, qb)
	# Previous step coordinates
	Δx1 = A * displacement(joint, xa1, qa1, xb1, qb1)

	# Finite difference
	Δv = (Δx - Δx1) / timestep
	return Δv
end
