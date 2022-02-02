################################################################################
# Test Derivatives
################################################################################

using Test

mech = get_pendulum()
eqc0 = mech.eqconstraints[1]
xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
λ = rand(2)
η = rand(0)
rot0 = eqc0.constraints[2]
GaT(rot0, xa, qa, xb, qb, η)
GbT(rot0, xa, qa, xb, qb, η)

# ∂aGb
J0 = ∂aGb(rot0, xa, qa, xb, qb, λ)
attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
J1 = FiniteDiff.finite_difference_jacobian(
    z -> GbT(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
    [xa; vector(qa)]
    ) * attjac
@test norm(J0 - J1, Inf) < 1e-7

# ∂bGb
J0 = ∂bGb(rot0, xa, qa, xb, qb, λ)
attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
J1 = FiniteDiff.finite_difference_jacobian(
    z -> GbT(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
    [xb; vector(qb)]
    ) * attjac
@test norm(J0 - J1, Inf) < 1e-7

# ∂aGa
J0 = ∂aGa(rot0, xa, qa, xb, qb, λ)
attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
J1 = FiniteDiff.finite_difference_jacobian(
    z -> GaT(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
    [xa; vector(qa)]
    ) * attjac
@test norm(J0 - J1, Inf) < 1e-7

# ∂bGa
J0 = ∂bGa(rot0, xa, qa, xb, qb, λ)
attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
J1 = FiniteDiff.finite_difference_jacobian(
    z -> GaT(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
    [xb; vector(qb)]
    ) * attjac
@test norm(J0 - J1, Inf) < 1e-7
