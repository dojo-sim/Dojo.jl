qrotate(q1::UnitQuaternion,q2::UnitQuaternion) = q2 * q1 / q2
vrotate(v::Vector,q::UnitQuaternion) = imag(qrotate(pure_quaternion(v), q))
vrotate(v::StaticVector,q::UnitQuaternion) = q*v

rotation_matrix(q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)

# ∂(rotation_matrix(q)*p)/∂q
∂qrotation_matrix(q::UnitQuaternion, p::AbstractVector) =
 	∂qVRᵀmat(LVᵀmat(q) * p) + VRᵀmat(q) * ∂qLVᵀmat(p)

# ∂(rotation_matrix(inv(q))*p)/∂q
∂qrotation_matrix_inv(q::UnitQuaternion, p::AbstractVector) =
 	∂qrotation_matrix(inv(q), p) * Tmat()

∂vrotate∂p(p::AbstractVector, q::UnitQuaternion) = VRᵀmat(q) * LVᵀmat(q)
∂vrotate∂q(p::AbstractVector, q::UnitQuaternion) = VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))
