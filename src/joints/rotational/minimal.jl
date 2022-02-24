################################################################################
# Displacements
################################################################################
function displacement(joint::Rotational,
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion;
        vmat=true)

    q = inv(joint.qoffset) * inv(qa) * qb
    vmat ? (return Vmat(q)) : (return q)
end

function displacement_jacobian_configuration(relative::Symbol, joint::Rotational,
        xa::AbstractVector{T}, qa::UnitQuaternion,
        xb::AbstractVector{T}, qb::UnitQuaternion;
        attjac::Bool=true, vmat=true) where T
    X = szeros(T, 3, 3)
    if relative == :parent
		Q = Lᵀmat(joint.qoffset) * Rmat(qb) * Tmat()
		attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
		Q = Lᵀmat(joint.qoffset) * Lᵀmat(qa)
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
    ∂rv∂q = ∂rotation_vector∂q(q)

	return A * [X ∂rv∂q * Q]
end

################################################################################
# Set Coordinates
################################################################################
function set_minimal_coordinates!(joint::Rotational, 
	pnode::Node, cnode::Node, 
	timestep;
	Δθ::AbstractVector=szeros(control_dimension(joint)))
	
	qoffset = joint.qoffset
	qa = pnode.state.q2[1]
	Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
	Δq = axis_angle_to_quaternion(Aᵀ * Δθ)
	qb = qa * qoffset * Δq
	set_maximal_configuration!(cnode; x=cnode.state.x2[1], q = qb)
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

	qoffset = joint.qoffset
	A = nullspace_mask(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(qoffset) * inv(qa) * qb
	q1 = inv(qoffset) * inv(qa1) * qb1

	return A * rotation_vector(inv(q1) * q) ./ timestep
end

function minimal_velocities_jacobian_configuration(relative::Symbol, joint::Rotational{T},
	xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
	xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
	timestep) where T
	
	qoffset = joint.qoffset
	A = nullspace_mask(joint)
	nu = control_dimension(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(qoffset) * inv(qa) * qb
	q1 = inv(qoffset) * inv(qa1) * qb1

	# return A * rotation_vector(inv(q1) * q) ./ timestep

	X = szeros(T, nu, 3)

	if relative == :parent 
		Q = 1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Rmat(q) * Tmat() * Rmat(qb1) * Lmat(inv(qoffset)) * Tmat() * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
		Q += 1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Lmat(inv(q1)) * Rmat(qb) * Lmat(inv(qoffset)) * Tmat()
		Q *= LVᵀmat(qa)
	elseif relative == :child 
		Q = 1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(qoffset) * inv(qa1)) * Rmat(quaternion_map(-ϕb, timestep)) * timestep / 2
		Q += 1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Lmat(inv(q1) * inv(qoffset) * inv(qa))
		Q *= LVᵀmat(qb)
	end

	return [X Q]
end

function minimal_velocities_jacobian_velocity(relative::Symbol, joint::Rotational{T},
	xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
	xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
	timestep) where T
	
	qoffset = joint.qoffset
	A = nullspace_mask(joint)
	nu = control_dimension(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	q = inv(qoffset) * inv(qa) * qb
	q1 = inv(qoffset) * inv(qa1) * qb1

	# return A * rotation_vector(inv(q1) * q) ./ timestep

	V = szeros(T, nu, 3)

	if relative == :parent
		Ω = -1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(qoffset)) * Rmat(qb1) * Tmat() * Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2 
	elseif relative == :child 
		Ω = -1.0 / timestep * A * ∂rotation_vector∂q(inv(q1) * q) * Rmat(q) * Tmat() * Lmat(inv(qoffset) * inv(qa1)) * Lmat(qb) * quaternion_map_jacobian(-ϕb, timestep) * timestep / 2
	end

	return [V Ω]
end
