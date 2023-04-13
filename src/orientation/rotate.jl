# rotate quaternion
quaternion_rotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2

# rotate vector
vector_rotate(v::AbstractVector,q::Quaternion) = Vmat(quaternion_rotate(Quaternion(0,v...), q))
∂vector_rotate∂q(p::AbstractVector, q::Quaternion) = VLmat(q) * Lmat(Quaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(Quaternion(p))

# rotate matrix
function matrix_rotate(A::AbstractMatrix,q::Quaternion)
    c1 = vector_rotate(A[SA[1;2;3]], q)
    c2 = vector_rotate(A[SA[4;5;6]], q)
    c3 = vector_rotate(A[SA[7;8;9]], q)

    return [c1 c2 c3]
end

# transform matrix in different frame "q*A*inv(q)"
function matrix_transform(A::AbstractMatrix,q::Quaternion)
    matrix_rotate(matrix_rotate(A', q)', q)
end

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
