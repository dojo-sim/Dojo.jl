using Dojo
using Dojo: params, pure_quaternion, Lmat, Lᵀmat, Rmat, Rᵀmat, Tmat, Tᵀmat, Vmat, Vᵀmat, VLmat, VLᵀmat, VRmat, VRᵀmat, LVᵀmat, LᵀVᵀmat, RVᵀmat, RᵀVᵀmat, slerp
using StaticArrays
using LinearAlgebra

w1 = rand()
v1 = rand(3)
w2 = rand()
v2 = rand(3)
qvref1 = pure_quaternion(v1)
qref1 = UnitQuaternion(w1,v1...)
qref2 = UnitQuaternion(w2,v2...)
q1 = UnitQuaternion(w1,SA[v1...])
@test q1 == qref1
q1 = UnitQuaternion(w1,v1)
@test q1 == qref1
qv1 = UnitQuaternion(SA[v1...])
@test qv1 == qvref1
qv1 = UnitQuaternion(v1)
@test qv1 == qvref1


@test isapprox(norm(params(qref1*qref2) - Lmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qref2) - Lᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2*qref1) - Rmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2/qref1) - Rᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)

@test isapprox(norm(params(inv(qref1)) - Tmat()*params(q1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(inv(qref1)) - Tᵀmat()*params(q1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1)[2:4] - Vmat()*params(qv1)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1) - Vᵀmat()*v1), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1*qref2)[2:4] - VLmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qref2)[2:4] - VLᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2*qref1)[2:4] - VRmat(q1)*params(qref2)), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2/qref1)[2:4] - VRᵀmat(q1)*params(qref2)), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1*qvref1) - LVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref1\qvref1) - LᵀVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1*qref1) - RVᵀmat(q1)*v1), 0.0; atol = 1e-10)
@test isapprox(norm(params(qvref1/qref1) - RᵀVᵀmat(q1)*v1), 0.0; atol = 1e-10)

@test isapprox(norm(params(qref1) - params(slerp(qref1,qref2,0.0))), 0.0; atol = 1e-10)
@test isapprox(norm(params(qref2) - params(slerp(qref1,qref2,1.0))), 0.0; atol = 1e-10)