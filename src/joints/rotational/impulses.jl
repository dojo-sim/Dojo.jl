################################################################################
 # Impulses
################################################################################
# see joints/joint.jl method

function impulse_map(relative::Symbol, joint::Rotational{T,Nλ,0},
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion,
    η) where {T,Nλ}
	return impulse_transform(relative, joint, xa, qa, xb, qb) * impulse_projector(joint)
end


################################################################################
 # Derivatives
################################################################################
function impulse_transform_jacobian(relative::Symbol, jacobian::Symbol,
    joint::Rotational{T,Nλ},
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion, p; attjac=true) where {T,Nλ}

    if relative == :parent
        if jacobian == :parent
            # ∂(Ja'*p)/∂(xa,qa)
            Z3 = szeros(T,3,3)
            # ∂Q∂qa = ∂VLᵀmat∂q(Tmat() * Rᵀmat(qb) * LVᵀmat(joint.axis_offset) * p) * LVᵀmat(qa)
            # return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∂Q∂qa]
            # ∂Q∂qa = Z3

            return szeros(T,6,6)
			# ∂Q∂qa = -∂rotation_matrix_inv∂q(qa, p, attjac=true)
			# return [Z3 Z3; Z3 ∂Q∂qa]
        elseif jacobian == :child
            # ∂(Ja'*p)/∂(xb,qb)
            Z3 = szeros(T,3,3)
            # ∂Q∂qb = VLᵀmat(qa) * Tmat(T) * ∂Rᵀmat∂q(LVᵀmat(joint.axis_offset) * p) * LVᵀmat(qb)
            # return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∂Q∂qb]
            # ∂Q∂qb = Z3

			return szeros(T,6,6)
        end
    elseif relative == :child
        if jacobian == :parent
            # ∂(Jb'*p)/∂(xa,qa)
            Z3 = szeros(T,3,3)
            # ∂Q∂qa = VLᵀmat(qb) * ∂Lmat∂q(LVᵀmat(joint.axis_offset) * p) * LVᵀmat(qa)
            # return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∂Q∂qa]

			∂Q∂qa = rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, p, attjac=true)
            # return [Z3 Z3; Z3 ∂Q∂qa]
			# return szeros(T,6,6)
        elseif jacobian == :child
        # ∂(Jb'*p)/∂(xb,qb)
            Z3 = szeros(T,3,3)
            # ∂Q∂qb = ∂VLᵀmat∂q(Lmat(qa) * LVᵀmat(joint.axis_offset) * p) * LVᵀmat(qb)
            # return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∂Q∂qb]

			∂Q∂qb = ∂rotation_matrix_inv∂q(qb, rotation_matrix(qa) * p, attjac=true)
            # ∂Q∂qb = ∂rotation_matrix_inv∂q(qb, p, attjac=true)
            return [Z3 Z3; Z3 ∂Q∂qb]
        end
    end
end
