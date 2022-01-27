################################################################################
# Test Derivatives
################################################################################
@testset "translational joint: impulse_map" begin
    mech = getpendulum()
    joint0 = mech.joints[1]
    tra0 = joint0.constraints[1]

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)
    λ = rand(3)
    η = rand(0)
    impulse_map_parent(tra0, xa, qa, xb, qb, η)
    impulse_map_child(tra0, xa, qa, xb, qb, η)

    # impulse_map_parent_jacobian_parent
    J0 = impulse_map_parent_jacobian_parent(tra0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_parent(tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_parent_jacobian_child
    J0 = impulse_map_parent_jacobian_child(tra0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_parent(tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_parent
    J0 = impulse_map_child_jacobian_parent(tra0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_child(tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_child
    J0 = impulse_map_child_jacobian_child(tra0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_child(tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end

@testset "rotational joint: impulse_map" begin
    mech = getpendulum()
    joint0 = mech.joints[1]
    rot0 = joint0.constraints[2]

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)
    λ = rand(2)
    η = rand(0)
    impulse_map_parent(rot0, xa, qa, xb, qb, η)
    impulse_map_child(rot0, xa, qa, xb, qb, η)

    # impulse_map_parent_jacobian_parent
    J0 = impulse_map_parent_jacobian_parent(rot0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_parent(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_parent_jacobian_child
    J0 = impulse_map_parent_jacobian_child(rot0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_parent(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_parent
    J0 = impulse_map_child_jacobian_parent(rot0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_child(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_child
    J0 = impulse_map_child_jacobian_child(rot0, xa, qa, xb, qb, λ)
    attjac = cat(I(3),LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> impulse_map_child(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end
