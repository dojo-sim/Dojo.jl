################################################################################
# Test Constraints
################################################################################

mech = get_snake(spring=10.0, damper=1.0, Nb=2, gravity=0.0, contact=false, jointtype=:Revolute)
initialize!(mech, :snake)
joint0 = mech.joints[2]
tra0 = joint0.constraints[1]
rot0 = joint0.constraints[2]

xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
η0 = 0
constraint(rot0, xa, qa, xb, qb, η0)

∇0 = constraint_jacobian_parent(rot0, xa, qa, xb, qb, η0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> constraint(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb, η0), [xa; vector(qa)])
     # * cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = constraint_jacobian_child(rot0, xa, qa, xb, qb, η0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> constraint(rot0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false), η0), [xb; vector(qb)])
     # * cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7


p0 = rand(3)
impulse_transform_parent(rot0, xa, qa, xb, qb)
impulse_transform_child(rot0, xa, qa, xb, qb)

∇0 = impulse_transform_parent_jacobian_parent(rot0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_parent(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb)*p0, [xa; vector(qa)]) *
    cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_parent_jacobian_child(rot0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_parent(rot0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false))*p0, [xb; vector(qb)]) *
    cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_child_jacobian_parent(rot0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_child(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb)*p0, [xa; vector(qa)]) *
    cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_child_jacobian_child(rot0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_child(rot0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false))*p0, [xb; vector(qb)]) *
    cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7



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
