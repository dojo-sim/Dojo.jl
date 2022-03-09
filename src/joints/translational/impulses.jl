################################################################################
 # Impulses
################################################################################
# see joints/joint.jl method 

################################################################################
# Impulse Transform Derivatives
################################################################################
function impulse_transform_jacobian(relative::Symbol, jacobian::Symbol,
    joint::Translational{T,Nλ},
    xa::AbstractVector, qa::Quaternion, 
    xb::AbstractVector, qb::Quaternion, p; 
    attjac=true) where {T,Nλ}

    Z3 = szeros(T,3,3)

    if relative == :parent 
        if jacobian == :parent 
            # ∂(impulse_transform_a'*p)/∂(xa,qa)
            ∇Xqa = -∂rotation_matrix∂q(qa, p) * LVᵀmat(qa)
            ∇Qxa =  ∂skew∂p(p) * rotation_matrix(inv(qa))
            ∇Qqa = -∂skew∂p(p) * ∂rotation_matrix_inv∂q(qa, xb - xa + rotation_matrix(qb) * joint.vertices[2]) * LVᵀmat(qa)
            return [Z3 ∇Xqa; ∇Qxa ∇Qqa]
        elseif jacobian == :child 
            # ∂(impulse_transform_a'*p)/∂(xb,qb)
            ∇Qxb = -∂skew∂p(p) * rotation_matrix(inv(qa))
            ∇Qqb = -∂skew∂p(p) * rotation_matrix(inv(qa)) * ∂rotation_matrix∂q(qb, joint.vertices[2]) * LVᵀmat(qb)
            return [Z3 Z3; ∇Qxb ∇Qqb]
        end
    elseif relative == :child 
        if jacobian == :parent 
            # ∂(impulse_transform_b'*p)/∂(xa,qa)
            cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
            ∇Xqa = ∂rotation_matrix∂q(qa, p) * LVᵀmat(qa)
            ∇Qqa = rotation_matrix(inv(qb)) * skew(cbpb_w) * ∂rotation_matrix∂q(qa, p) * LVᵀmat(qa)
            return [Z3 ∇Xqa; Z3 ∇Qqa]
        elseif jacobian == :child 
            # ∂(impulse_transform_b'*p)/∂(xb,qb)
            cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
            ∇Qqb = ∂rotation_matrix_inv∂q(qb, skew(cbpb_w) * rotation_matrix(qa) * p)
            ∇Qqb += rotation_matrix(inv(qb)) * ∂skew∂p(rotation_matrix(qa) * p) * ∂rotation_matrix∂q(qb, joint.vertices[2])
            ∇Qqb *= LVᵀmat(qb)
            return [Z3 Z3; Z3 ∇Qqb]
        end
    end
end