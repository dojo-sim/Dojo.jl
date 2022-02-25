################################################################################
# Displacements
################################################################################
function displacement(joint::Rotational,
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion;
        vmat=true)

    q = inv(joint.axis_offset) * inv(qa) * qb
    vmat ? (return Vmat(q)) : (return q)
end

function displacement_jacobian_configuration(relative::Symbol, joint::Rotational,
        xa::AbstractVector{T}, qa::UnitQuaternion,
        xb::AbstractVector{T}, qb::UnitQuaternion;
        attjac::Bool=true, vmat=true) where T
    X = szeros(T, 3, 3)
    if relative == :parent
		Q = Lᵀmat(joint.axis_offset) * Rmat(qb) * Tmat()
		attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
		Q = Lᵀmat(joint.axis_offset) * Lᵀmat(qa)
		attjac && (Q *= LVᵀmat(qb))
	end
	vmat && (Q = Vmat() * Q)
	return X, Q
end

################################################################################
# Coordinates
################################################################################
function minimal_coordinates(joint::Rotational,
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion)

    return nullspace_mask(joint) * rotation_vector(displacement(joint, xa, qa, xb, qb, vmat=false))
end

function minimal_coordinates_jacobian_configuration(relative::Symbol, joint::Rotational{T},
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion;
        attjac::Bool=true) where T

    A = nullspace_mask(joint)
    q = displacement(joint, xa, qa, xb, qb, vmat=false)
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=attjac, vmat=false)
    ∂rv∂q = drotation_vectordq(q)

	return A * [X ∂rv∂q * Q]
end

################################################################################
# Set Coordinates
################################################################################
function set_minimal_coordinates!(joint::Rotational, 
	pnode::Node, cnode::Node, 
	timestep;
	Δθ::AbstractVector=szeros(input_dimension(joint)))
	
	axis_offset = joint.axis_offset
	qa = pnode.state.q2
	Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
	Δq = axis_angle_to_quaternion(Aᵀ * Δθ)
	qb = qa * axis_offset * Δq
	set_maximal_coordinates!(cnode; x=cnode.state.x2, q = qb)
	set_previous_configuration!(cnode, timestep)
	return nothing
end

################################################################################
# Velocities
################################################################################
function minimal_velocities(joint::Rotational,
		xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
		xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
		timestep)

	axis_offset = joint.axis_offset
	A = nullspace_mask(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(axis_offset) * inv(qa) * qb
	q1 = inv(axis_offset) * inv(qa1) * qb1

	return A * rotation_vector(inv(q1) * q) ./ timestep
end

function minimal_velocities_jacobian_configuration(relative::Symbol, joint::Rotational{T},
	xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
	xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
	timestep) where T
	
	axis_offset = joint.axis_offset
	A = nullspace_mask(joint)
	nu = input_dimension(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(axis_offset) * inv(qa) * qb
	q1 = inv(axis_offset) * inv(qa1) * qb1

	X = szeros(T, nu, 3)

	if relative == :parent 
		Q = 1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Rmat(q) * Tmat() * Rmat(qb1) * Lmat(inv(axis_offset)) * Tmat() * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
		Q += 1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Lmat(inv(q1)) * Rmat(qb) * Lmat(inv(axis_offset)) * Tmat()
		Q *= LVᵀmat(qa)
	elseif relative == :child 
		Q = 1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(axis_offset) * inv(qa1)) * Rmat(quaternion_map(-ϕb, timestep)) * timestep / 2
		Q += 1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Lmat(inv(q1) * inv(axis_offset) * inv(qa))
		Q *= LVᵀmat(qb)
	end

	return [X Q]
end

function minimal_velocities_jacobian_velocity(relative::Symbol, joint::Rotational{T},
	xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
	xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
	timestep) where T
	
	axis_offset = joint.axis_offset
	A = nullspace_mask(joint)
	nu = input_dimension(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(axis_offset) * inv(qa) * qb
	q1 = inv(axis_offset) * inv(qa1) * qb1

	V = szeros(T, nu, 3)

	if relative == :parent
		Ω = -1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(axis_offset)) * Rmat(qb1) * Tmat() * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
	elseif relative == :child 
		Ω = -1.0 / timestep * A * drotation_vectordq(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(axis_offset) * inv(qa1)) * Lmat(qb) * quaternion_map_jacobian(-ϕb, timestep) * timestep / 2
	end

	return [V Ω]
end
