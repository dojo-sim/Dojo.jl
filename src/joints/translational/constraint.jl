mutable struct Translational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N,Nb½}
    axis::SVector{3,T} # translation axis in parent frame
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    spring_type::Symbol # the rotational springs can be :sinusoidal or :linear, if linear then we need joint_limits to avoid the 180° singularity.
    input::SVector{3,T}
end

function Translational{T,Nλ}(body1::Node, body2::Node;
        p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3),
        spring = zero(T), damper = zero(T), spring_offset = szeros(T,3-Nλ),
        joint_limits = [szeros(T,0), szeros(T,0)],
        spring_type::Symbol = :sinusoidal,
    ) where {T,Nλ}
    vertices = (p1, p2)
    V1, V2, V3 = orthogonal_rows(axis)
    V12 = [V1;V2]
    input = zeros(T,3)
    Nb½ = length(joint_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    Translational{T,Nλ,Nb,N,Nb½,N̄λ}(axis, V3, V12, vertices, spring, damper, spring_offset, joint_limits, spring_type, input), body1.id, body2.id
end

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

# function impulse_transform_jacobian(:parent, :child, joint::Translational{T,Nλ},
#         xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
#     # ∂(impulse_transform_a'*p)/∂(xb,qb)
#     Z3 = szeros(T,3,3)

#     ∇Qxb = -∂pskew(p) * rotation_matrix(inv(qa))
#     ∇Qqb = -∂pskew(p) * rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, joint.vertices[2]) * LVᵀmat(qb)
#     return [Z3   Z3; ∇Qxb ∇Qqb]
# end

# function impulse_transform_jacobian(:child, :parent, joint::Translational{T,Nλ},
#         xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
#     # ∂(impulse_transform_b'*p)/∂(xa,qa)
#     Z3 = szeros(T,3,3)
#     cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point

#     ∇Xqa = ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
#     ∇Qqa = rotation_matrix(inv(qb)) * skew(cbpb_w) * ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
#     return [Z3 ∇Xqa; Z3 ∇Qqa]
# end

# function impulse_transform_jacobian(:child, :child, joint::Translational{T,Nλ},
#         xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
#     # ∂(impulse_transform_b'*p)/∂(xb,qb)
#     Z3 = szeros(T,3,3)
#     cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
#     ∇Qqb = ∂qrotation_matrix_inv(qb, skew(cbpb_w) * rotation_matrix(qa) * p)
#     ∇Qqb += rotation_matrix(inv(qb)) * ∂pskew(rotation_matrix(qa) * p) * ∂qrotation_matrix(qb, joint.vertices[2])
#     ∇Qqb *= LVᵀmat(qb)
#     return [Z3 Z3; Z3 ∇Qqb]
# end
