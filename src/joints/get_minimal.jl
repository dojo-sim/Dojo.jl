################################################################################
# GET: Coordinates Joint
################################################################################
@inline function minimal_coordinates(joint::Joint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::JointConstraint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
	Δx = minimal_coordinates(joint.translational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
    Δθ = minimal_coordinates(joint.rotational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
	return [Δx; Δθ]
end


################################################################################
# GET: Coordinates Joint
################################################################################
@inline function minimal_coordinates_velocities(joint::JointConstraint,
		pnode::Node, cnode::Node, timestep)
	Δxθ = minimal_coordinates(joint, pnode, cnode)
	Δvϕ = minimal_velocities(joint, pnode, cnode, timestep)
	return [Δxθ; Δvϕ]
end

@inline function minimal_velocities(joint::JointConstraint, pnode::Node, cnode::Node, timestep)
	Δv = minimal_velocities(joint.translational, pnode, cnode, timestep)
	Δϕ = minimal_velocities(joint.rotational, pnode, cnode, timestep)
	return [Δv; Δϕ]
end

@inline function minimal_velocities(joint::Joint, pnode::Node, cnode::Node, timestep)
	minimal_velocities(joint, initial_configuration_velocity(pnode.state)...,
		initial_configuration_velocity(cnode.state)..., timestep)
end

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

@inline function minimal_velocities(joint::Rotational,
		xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
		xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
		timestep)
	qoffset = joint.qoffset
	A = nullspace_mask(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	# Coordinates
	q = inv(qoffset) * inv(qa) * qb
	Δθ = A * rotation_vector(q)
	# Previous step coordinates
	q1 = inv(qoffset) * inv(qa1) * qb1
	Δθ1 = A * rotation_vector(q1)

	# Finite difference
    Δϕ = (Δθ - Δθ1) / timestep
	return Δϕ
end

@inline function minimal_velocities_jacobian_configuration(relative::Symbol, joint::Joint{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector, timestep) where T

    if relative == :parent
		∇xq = FiniteDiff.finite_difference_jacobian(xq -> minimal_velocities(
			joint, xq[SUnitRange(1,3)], va, UnitQuaternion(xq[4:7]..., false),
			ϕa, xb, vb, qb, ϕb, timestep), [xa; vector(qa)]) * cat(I(3), LVᵀmat(qa), dims=(1,2))
    elseif relative == :child
		∇xq = FiniteDiff.finite_difference_jacobian(xq -> minimal_velocities(
			joint, xa, va, qa, ϕa, xq[SUnitRange(1,3)], vb, UnitQuaternion(xq[4:7]..., false),
			ϕb, timestep), [xb; vector(qb)]) * cat(I(3), LVᵀmat(qb), dims=(1,2))
    end
    return ∇xq
end

@inline function minimal_velocities_jacobian_velocity(relative::Symbol, joint::Joint{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector, timestep) where T

	if relative == :parent
		∇vϕ = FiniteDiff.finite_difference_jacobian(vϕ -> minimal_velocities(
			joint, xa, vϕ[SUnitRange(1,3)], qa, vϕ[SUnitRange(4,6)],
			xb, vb, qb, ϕb, timestep), [va; ϕa])
    elseif relative == :child
		∇vϕ = FiniteDiff.finite_difference_jacobian(vϕ -> minimal_velocities(
			joint, xa, va, qa, ϕa, xb, vϕ[SUnitRange(1,3)], qb, vϕ[SUnitRange(4,6)],
			timestep), [vb; ϕb])
    end
    return ∇vϕ
end
