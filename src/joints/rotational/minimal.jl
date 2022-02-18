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
    # ∂rv∂q = FiniteDiff.finite_difference_jacobian(r -> rotation_vector(UnitQuaternion(r..., false)), vector(q))
    ∂rv∂q = ∂rotation_vector∂q(q)

	return A * [X ∂rv∂q * Q]
end

################################################################################
# Set Coordinates
################################################################################
function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::Rotational, timestep;
	Δθ::AbstractVector=szeros(control_dimension(joint)))
	# Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
	qoffset = joint.qoffset
	qa = pnode.state.q2[1]
	Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
	Δq = axis_angle_to_quaternion(Aᵀ*Δθ)
	qb = qa * qoffset * Δq
	set_position!(cnode; x=cnode.state.x2[1], q = qb)
	set_previous_configuration!(cnode, timestep)
	return nothing
end



################################################################################
# Velocities
################################################################################
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

	# # Coordinates
	# q = inv(qoffset) * inv(qa) * qb
	# Δθ = A * rotation_vector(q)
	# # Previous step coordinates
	# q1 = inv(qoffset) * inv(qa1) * qb1
	# Δθ1 = A * rotation_vector(q1)
	#
	# # Finite difference
    # Δϕ = (Δθ - Δθ1) / timestep



	q = inv(qoffset) * inv(qa) * qb
    q1 = inv(qoffset) * inv(qa1) * qb1
    @show sol = A * rotation_vector(inv(q1) * q) ./ timestep
    return sol

	# return Δϕ
end
