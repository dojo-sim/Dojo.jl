# rotate quaternion
quaternion_rotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2

# rotate vector
vector_rotate(v::AbstractVector,q::Quaternion) = Vmat(quaternion_rotate(Quaternion(v), q))
# vector_rotate(v::StaticVector,q::Quaternion) = q * v
∂vector_rotate∂q(p::AbstractVector, q::Quaternion) = VLmat(q) * Lmat(Quaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(Quaternion(p))

# rotation matrix
rotation_matrix(q::Quaternion) = VRᵀmat(q) * LVᵀmat(q)
function ∂rotation_matrix∂q(q::Quaternion, p::AbstractVector{T}) where T #; attjac::Bool=false) where T
    M = ∂VRᵀmat∂q(LVᵀmat(q) * p) + VRᵀmat(q) * ∂LVᵀmat∂q(p)
    # if attjac
    #     return (M * LVᵀmat(q))::SMatrix{3,3,T,9}
    # else
    #     return M::SMatrix{3,4,T,12}
    # end
end
function ∂rotation_matrix_inv∂q(q::Quaternion, p::AbstractVector{T}) where T #; attjac::Bool=false) where T
    M = ∂rotation_matrix∂q(inv(q), p) * Tmat()
    # if attjac
    #     return (M * LVᵀmat(q))::SMatrix{3,3,T,9}
    # else
    #     return M::SMatrix{3,4,T,12}
    # end
end
