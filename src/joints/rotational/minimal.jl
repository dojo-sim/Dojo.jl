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
    ∂rotation_vector∂q = FiniteDiff.finite_difference_jacobian(r -> rotation_vector(UnitQuaternion(r..., false)), vector(q))

	return A * [X ∂rotation_vector∂q * Q]
end
