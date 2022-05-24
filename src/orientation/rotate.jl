# rotate quaternion
quaternion_rotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2

# rotate vector
vector_rotate(v::AbstractVector,q::Quaternion) = Vmat(quaternion_rotate(Quaternion(v), q))
∂vector_rotate∂q(p::AbstractVector, q::Quaternion) = VLmat(q) * Lmat(Quaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(Quaternion(p))

# # rotate matrix
# function matrix_rotate(A::AbstractMatrix,q::Quaternion)
#     c1 = Vmat(quaternion_rotate(Quaternion(A[SVector{3}(1,2,3)]), q))
#     c2 = Vmat(quaternion_rotate(Quaternion(A[SVector{3}(4,5,6)]), q))
#     c3 = Vmat(quaternion_rotate(Quaternion(A[SVector{3}(7,8,9)]), q))

#     return svcat(c1, c2, c3)
# end

# rotation matrix
rotation_matrix(q::Quaternion) = VRᵀmat(q) * LVᵀmat(q)
function ∂rotation_matrix∂q(q::Quaternion, p::AbstractVector; attjac::Bool=false)
    M = ∂VRᵀmat∂q(LVᵀmat(q) * p) + VRᵀmat(q) * ∂LVᵀmat∂q(p)
    attjac && (M *= LVᵀmat(q))
    return M
end
function ∂rotation_matrix_inv∂q(q::Quaternion, p::AbstractVector; attjac::Bool=false)
    M = ∂rotation_matrix∂q(inv(q), p) * Tmat()
    attjac && (M *= LVᵀmat(q))
    return M
end
