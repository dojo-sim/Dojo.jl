################################################################################
 # Impulses
################################################################################
# see joints/joint.jl method 

################################################################################
 # Derivatives
################################################################################
function impulse_transform_jacobian(relative::Symbol, jacobian::Symbol,
    joint::Rotational{T,Nλ},
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, p; attjac=true) where {T,Nλ}

    if relative == :parent 
        if jacobian == :parent 
            # ∂(Ja'*p)/∂(xa,qa)
            Z3 = szeros(T,3,3)
            ∂Q∂qa = ∂qVLᵀmat(Tmat() * Rᵀmat(qb) * LVᵀmat(joint.qoffset) * p) * LVᵀmat(qa)
            return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∂Q∂qa]
        elseif jacobian == :child 
            # ∂(Ja'*p)/∂(xb,qb)
            Z3 = szeros(T,3,3)
            ∇Qqb = VLᵀmat(qa) * Tmat(T) * ∂qRᵀmat(LVᵀmat(joint.qoffset) * p) * LVᵀmat(qb)
            return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∇Qqb]
        end
    elseif relative == :child 
        if jacobian == :parent 
            # ∂(Jb'*p)/∂(xa,qa)
            Z3 = szeros(T,3,3)
            ∇Qqa = VLᵀmat(qb) * ∂qLmat(LVᵀmat(joint.qoffset) * p) * LVᵀmat(qa)
            return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∇Qqa]
        elseif jacobian == :child 
        # ∂(Jb'*p)/∂(xb,qb)
            Z3 = szeros(T,3,3)
            ∇Qqb = ∂qVLᵀmat(Lmat(qa) * LVᵀmat(joint.qoffset) * p) * LVᵀmat(qb)
            return cat(I(3), 0.5 * I(3), dims=(1,2)) * [Z3 Z3; Z3 ∇Qqb]
        end
    end
end