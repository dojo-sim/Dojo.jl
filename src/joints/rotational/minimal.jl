@inline function orientation_error(joint::Rotational, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    # (norm(Vmat(joint.qoffset)) > 1e-6) && (@warn "check the validity of this expression")

    q = inv(joint.qoffset) * inv(qa) * qb # rotation from frame b to frame a then joint_qoffset then spring_qoffset
    "q = Lᵀmat(qoffset) * Rmat(qb) * Tmat() * vector(qa)"
    "q = Lᵀmat(qoffset) * Lᵀmat(qa) * vector(qb)"
    # b --qb-> world --inv(qa)-> a --inv(qoffset)-> jointoffsetframe --inv(qoffset)-> springoffsetframe
    return q
end

function orientation_error_jacobian_configuration(jacobian_relative::Symbol,
        joint::Rotational, xa::AbstractVector{T}, qa::UnitQuaternion,
        xb::AbstractVector{T}, qb::UnitQuaternion; attjac::Bool=true) where T
	
    X = szeros(T, 3, 3)

    if jacobian_relative == :parent
		Q = Lᵀmat(joint.qoffset) * Rmat(qb) * Tmat()
		attjac && (Q *= LVᵀmat(qa))
    elseif jacobian_relative == :child
		Q = Lᵀmat(joint.qoffset) * Lᵀmat(qa)
		attjac && (Q *= LVᵀmat(qb))
	end

    return X, Q
end

################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Rotational, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    q = inv(joint.qoffset) * inv(qa) * qb
    return nullspace_mask(joint) * rotation_vector(q)
end

@inline function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol,
        joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector,
        qb::UnitQuaternion; attjac::Bool=true) where T
    A = nullspace_mask(joint)
    q = inv(joint.qoffset) * inv(qa) * qb
    X = szeros(T,3,3)
    if jacobian_relative == :parent
		# Q = ∂qrotation_vector(q) * Lᵀmat(joint.qoffset) * Rmat(qb) * Tmat()
		Q = FiniteDiff.finite_difference_jacobian(qa -> rotation_vector(inv(joint.qoffset) * inv(UnitQuaternion(qa..., false)) * qb), vector(qa))
		attjac && (Q *= LVᵀmat(qa))
    elseif jacobian_relative == :child
		# Q = ∂qrotation_vector(q) * Lᵀmat(joint.qoffset) * Lᵀmat(qa)
		# Q = ∂qrotation_vector(q) * rotation_matrix(inv(joint.qoffset) * inv(qa)) * ∂qrotation_matrix(qb, )
		Q = FiniteDiff.finite_difference_jacobian(qb -> rotation_vector(inv(joint.qoffset) * inv(qa) * UnitQuaternion(qb..., false)), vector(qb))
		attjac && (Q *= LVᵀmat(qb))
    end
	return A * [X Q]
end


################################################################################
# Velocities
################################################################################
@inline function minimal_velocities(joint::Rotational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
	qoffset = joint.qoffset
    Δϕ = rotation_matrix(inv(qoffset)) * (vrotate(ϕb, inv(qa) * qb) - ϕa) # in offset frame
    return nullspace_mask(joint) * Δϕ
end

@inline function minimal_velocities_jacobian_configuration_parent(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    X = szeros(T,3,3)
    Q = rotation_matrix(inv(qoffset)) * ∂qrotation_matrix_inv(qa, vrotate(ϕb, qb))
    ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qa)]
    return ∇xq
end
@inline function minimal_velocities_jacobian_configuration_child(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    X = szeros(T,3,3)
    Q = rotation_matrix(inv(qoffset) * inv(qa)) * ∂qrotation_matrix(qb, ϕb)
    ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qb)]
    return ∇xq
end
@inline function minimal_velocities_jacobian_velocity_parent(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    V = szeros(T,3,3)
    Ω = - rotation_matrix(inv(qoffset))
    ∇vϕ = nullspace_mask(joint) * [V Ω]
    return ∇vϕ
end
@inline function minimal_velocities_jacobian_velocity_child(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    V = szeros(T,3,3)
    Ω = rotation_matrix(inv(qoffset) * inv(qa) * qb)
    ∇vϕ = nullspace_mask(joint) * [V Ω]
    return ∇vϕ
end
