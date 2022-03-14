# rotate quaternion
quaternion_rotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2

# rotate vector
vector_rotate(v::Vector,q::UnitQuaternion) = imag(quaternion_rotate(pure_quaternion(v), q))
vector_rotate(v::StaticVector,q::UnitQuaternion) = q * v
∂vector_rotate∂p(p::AbstractVector, q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
∂vector_rotate∂q(p::AbstractVector, q::UnitQuaternion) = VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))

# rotation matrix
rotation_matrix(q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
∂rotation_matrix∂q(q::UnitQuaternion, p::AbstractVector) = ∂VRᵀmat∂q(LVᵀmat(q) * p) + VRᵀmat(q) * ∂LVᵀmat∂q(p)
∂rotation_matrix_inv∂q(q::UnitQuaternion, p::AbstractVector) = ∂rotation_matrix∂q(inv(q), p) * Tmat()

