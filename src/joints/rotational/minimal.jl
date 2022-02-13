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

################################################################################
# Velocities
################################################################################
function minimal_velocities(joint::Rotational,
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)

	qoffset = joint.qoffset
    Δϕ = rotation_matrix(inv(qoffset)) * (vrotate(ϕb, inv(qa) * qb) - ϕa) # in offset frame
    return nullspace_mask(joint) * Δϕ
end

function minimal_velocities_jacobian_configuration(relative::Symbol, joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T

    if relative == :parent
        qoffset = joint.qoffset
        X = szeros(T,3,3)
        Q = rotation_matrix(inv(qoffset)) * ∂qrotation_matrix_inv(qa, vrotate(ϕb, qb))
        ∇xq = nullspace_mask(joint) * [X Q * LVᵀmat(qa)]
    elseif relative == :child
        qoffset = joint.qoffset
        X = szeros(T,3,3)
        Q = rotation_matrix(inv(qoffset) * inv(qa)) * ∂qrotation_matrix(qb, ϕb)
        ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qb)]
    end

    return ∇xq
end

function minimal_velocities_jacobian_velocity(relative::Symbol, joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T

    if relative == :parent
        qoffset = joint.qoffset
        V = szeros(T,3,3)
        Ω = - rotation_matrix(inv(qoffset))
        ∇vϕ = nullspace_mask(joint) * [V Ω]
    elseif relative == :child
        qoffset = joint.qoffset
        V = szeros(T,3,3)
        Ω = rotation_matrix(inv(qoffset) * inv(qa) * qb)
        ∇vϕ = nullspace_mask(joint) * [V Ω]
    end

    return ∇vϕ
end
