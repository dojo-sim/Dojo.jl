# rotate quaternion
quaternion_rotate(q1::Quaternion,q2::Quaternion) = q2 * q1 / q2

# rotate vector
vector_rotate(v::AbstractVector,q::Quaternion) = Vmat(quaternion_rotate(Quaternion(v), q))
# vector_rotate(v::StaticVector,q::Quaternion) = q*v
∂vector_rotate∂p(p::AbstractVector, q::Quaternion) = VRᵀmat(q) * LVᵀmat(q)
∂vector_rotate∂q(p::AbstractVector, q::Quaternion) = VLmat(q) * Lmat(Quaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(Quaternion(p))

# rotation matrix
rotation_matrix(q::Quaternion) = VRᵀmat(q) * LVᵀmat(q)
∂rotation_matrix∂q(q::Quaternion, p::AbstractVector) = ∂VRᵀmat∂q(LVᵀmat(q) * p) + VRᵀmat(q) * ∂LVᵀmat∂q(p)
∂rotation_matrix_inv∂q(q::Quaternion, p::AbstractVector) = ∂rotation_matrix∂q(inv(q), p) * Tmat()

