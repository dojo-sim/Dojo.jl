################################################################################
# Test Constraints
################################################################################

mech = get_snake(spring=10.0, damper=1.0, Nb=2, gravity=0.0, contact=false, jointtype=:Revolute)
initialize!(mech, :snake)
mech = get_atlas()
initialize_atlas!(mech)
joint0 = mech.joints[2]
tra0 = joint0.constraints[1]
rot0 = joint0.constraints[2]

## rotation
xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
η0 = 0
constraint(rot0, xa, qa, xb, qb, η0)

∇0 = constraint_jacobian_configuration(:parent, rot0, xa, qa, xb, qb, η0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> constraint(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb, η0), [xa; vector(qa)])
     # * cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = constraint_jacobian_configuration(:child, rot0, xa, qa, xb, qb, η0)
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

## translation
xa = rand(3)
qa = UnitQuaternion(rand(4)...)
xb = rand(3)
qb = UnitQuaternion(rand(4)...)
η0 = 0
constraint(tra0, xa, qa, xb, qb, η0)

∇0 = constraint_jacobian_configuration(:parent, tra0, xa, qa, xb, qb, η0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> constraint(tra0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb, η0), [xa; vector(qa)])
     # * cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = constraint_jacobian_child(tra0, xa, qa, xb, qb, η0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> constraint(tra0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false), η0), [xb; vector(qb)])
     # * cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7


p0 = rand(3)
impulse_transform_parent(tra0, xa, qa, xb, qb)
impulse_transform_parent_alt(tra0, xa, qa, xb, qb)

impulse_transform_child(tra0, xa, qa, xb, qb)
0.5 * impulse_transform_child_alt(tra0, xa, qa, xb, qb)

∇0 = impulse_transform_parent_jacobian_parent(tra0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_parent(tra0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb)*p0, [xa; vector(qa)]) *
    cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_parent_jacobian_child(tra0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_parent(tra0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false))*p0, [xb; vector(qb)]) *
    cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_child_jacobian_parent(tra0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_child(tra0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb)*p0, [xa; vector(qa)]) *
    cat(I(3), LVᵀmat(qa), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7

∇0 = impulse_transform_child_jacobian_child(tra0, xa, qa, xb, qb, p0)
∇1 = FiniteDiff.finite_difference_jacobian(
    xq -> impulse_transform_child(tra0, xa, qa, xq[1:3], UnitQuaternion(xq[4:7]..., false))*p0, [xb; vector(qb)]) *
    cat(I(3), LVᵀmat(qb), dims=(1,2))
@test norm(∇0 - ∇1, Inf) < 1e-7


displacement(tra0, xa, qa, xb, qb, rotate=true)
displacement_jacobian_configuration(:parent, tra0, xa, qa, xb, qb)
Xa = FiniteDiff.finite_difference_jacobian(a -> displacement(tra0, a, qa, xb, qb, rotate=true), xa)
Qa = FiniteDiff.finite_difference_jacobian(a -> displacement(tra0, xa, UnitQuaternion(a..., false), xb, qb, rotate=true), vector(qa)) * LVᵀmat(qa)
transpose(impulse_transform_parent(tra0, xa, qa, xb, qb))

-VLᵀmat(qa) * RVᵀmat(qa)
(2.0 * VLᵀmat(qa) * (Lmat(UnitQuaternion(xb + vrotate(joint0.constraints[1].vertices[2], qb))) - Lmat(UnitQuaternion(xa)))) * LVᵀmat(qa)

displacement_jacobian_configuration(:child, tra0, xa, qa, xb, qb)
Xb = FiniteDiff.finite_difference_jacobian(b -> displacement(tra0, xa, qa, b, qb, rotate=true), xb)
Qb = FiniteDiff.finite_difference_jacobian(b -> displacement(tra0, xa, qa, xb, UnitQuaternion(b..., false), rotate=true), vector(qb)) * LVᵀmat(qb)
transpose(impulse_transform_child(tra0, xa, qa, xb, qb))

