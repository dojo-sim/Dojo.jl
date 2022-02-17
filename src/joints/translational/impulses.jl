################################################################################
 # Impulses
################################################################################
# see joints/joint.jl method 

################################################################################
# Impulse Transform Derivatives
################################################################################
function impulse_transform_jacobian(relative::Symbol, jacobian::Symbol,
    joint::Translational{T,Nλ},
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}

    Z3 = szeros(T,3,3)

    if relative == :parent 
        if jacobian == :parent 
            # ∂(impulse_transform_a'*p)/∂(xa,qa)
            ∇Xqa = -∂qrotation_matrix(qa, p) * LVᵀmat(qa)
            ∇Qxa =  ∂pskew(p) * rotation_matrix(inv(qa))
            ∇Qqa = -∂pskew(p) * ∂qrotation_matrix_inv(qa, xb - xa + rotation_matrix(qb) * joint.vertices[2]) * LVᵀmat(qa)
            return [Z3 ∇Xqa; ∇Qxa ∇Qqa]
        elseif jacobian == :child 
            # ∂(impulse_transform_a'*p)/∂(xb,qb)
            ∇Qxb = -∂pskew(p) * rotation_matrix(inv(qa))
            ∇Qqb = -∂pskew(p) * rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, joint.vertices[2]) * LVᵀmat(qb)
            return [Z3   Z3; ∇Qxb ∇Qqb]
        end
    elseif relative == :child 
        if jacobian == :parent 
            # ∂(impulse_transform_b'*p)/∂(xa,qa)
            cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
            ∇Xqa = ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
            ∇Qqa = rotation_matrix(inv(qb)) * skew(cbpb_w) * ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
            return [Z3 ∇Xqa; Z3 ∇Qqa]
        elseif jacobian == :child 
            # ∂(impulse_transform_b'*p)/∂(xb,qb)
            cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
            ∇Qqb = ∂qrotation_matrix_inv(qb, skew(cbpb_w) * rotation_matrix(qa) * p)
            ∇Qqb += rotation_matrix(inv(qb)) * ∂pskew(rotation_matrix(qa) * p) * ∂qrotation_matrix(qb, joint.vertices[2])
            ∇Qqb *= LVᵀmat(qb)
            return [Z3 Z3; Z3 ∇Qqb]
        end
    end
end