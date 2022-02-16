################################################################################
# Displacement Jacobian
################################################################################
@testset "displacement jacobian: rotational" begin
    mech = get_mechanism(:pendulum)
    joint0 = mech.joints[1]
    rot0 = joint0.rotational
    rot0.qoffset = UnitQuaternion(rand(4)...)


    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)

    # displacement_jacobian_configuration
    X0, Q0 = displacement_jacobian_configuration(:parent, rot0, xa, qa, xb, qb, attjac=true)
    J0 = [X0 Q0]
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.displacement(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb),
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # displacement_jacobian_configuration
    X0, Q0 = displacement_jacobian_configuration(:child, rot0, xa, qa, xb, qb, attjac=true)
    J0 = [X0 Q0]
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.displacement(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false)),
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end

@testset "displacement jacobian: translational" begin
    mech = get_mechanism(:slider)
    joint0 = mech.joints[1]
    tra0 = joint0.translational

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)

    # displacement_jacobian_configuration
    X0, Q0 = displacement_jacobian_configuration(:parent, tra0, xa, qa, xb, qb, attjac=true)
    J0 = [X0 Q0]
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.displacement(tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb),
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # displacement_jacobian_configuration
    X0, Q0 = displacement_jacobian_configuration(:child, tra0, xa, qa, xb, qb, attjac=true)
    J0 = [X0 Q0]
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.displacement(tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false)),
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end


################################################################################
# Impulse transform Jacobian
################################################################################
@testset "impulse transform jacobian: rotational" begin
    mech = get_mechanism(:pendulum)
    joint0 = mech.joints[1]
    rot0 = joint0.rotational
    rot0.qoffset = UnitQuaternion(rand(4)...)

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)

    p0 = rand(3)
    J0 = impulse_transform_parent_jacobian_parent(rot0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_parent(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb)*p0,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_parent_jacobian_child(rot0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_parent(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false))*p0,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_child_jacobian_parent(rot0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_child(rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb)*p0,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_child_jacobian_child(rot0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_child(rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false))*p0,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end

@testset "impulse transform jacobian: translational" begin
    mech = get_mechanism(:slider)
    joint0 = mech.joints[1]
    tra0 = joint0.translational

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)

    p0 = rand(3)
    J0 = impulse_transform_parent_jacobian_parent(tra0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_parent(tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb)*p0,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_parent_jacobian_child(tra0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_parent(tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false))*p0,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_child_jacobian_parent(tra0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_child(tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb)*p0,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    J0 = impulse_transform_child_jacobian_child(tra0, xa, qa, xb, qb, p0)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_transform_child(tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false))*p0,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end


################################################################################
# Impulse map Jacobian
################################################################################
@testset "impulse map: translational" begin
    mech = get_mechanism(:pendulum)
    joint0 = mech.joints[1]
    tra0 = joint0.translational

    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)
    λ = rand(3)
    η = rand(0)
    Dojo.impulse_map(:parent, tra0, xa, qa, xb, qb, η)
    Dojo.impulse_map(:child, tra0, xa, qa, xb, qb, η)

    m = 1.0
    J = I(3)
    bodya = Body(m, J)
    bodyb = Body(m, J)
    bodya.state.x2[1] = xa
    bodya.state.q2[1] = qa
    bodyb.state.x2[1] = xb
    bodyb.state.q2[1] = qb
    # impulse_map_parent_jacobian_parent
    J0 = Dojo.impulse_map_jacobian(:parent, :parent, tra0, bodya, bodyb, λ)

    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:parent, tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_parent_jacobian_child
    J0 = Dojo.impulse_map_jacobian(:parent, :child, tra0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:parent, tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_parent
    J0 = Dojo.impulse_map_jacobian(:child, :parent, tra0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:child, tra0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_child
    J0 = Dojo.impulse_map_jacobian(:child, :child, tra0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:child, tra0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end

@testset "impulse map: rotational" begin
    mech = get_mechanism(:pendulum)
    joint0 = mech.joints[1]
    rot0 = joint0.rotational
    rot0.qoffset = UnitQuaternion(rand(4)...)


    xa = rand(3)
    qa = UnitQuaternion(rand(4)...)
    xb = rand(3)
    qb = UnitQuaternion(rand(4)...)
    λ = rand(2)
    η = rand(0)
    Dojo.impulse_map(:parent, rot0, xa, qa, xb, qb, η)
    Dojo.impulse_map(:child, rot0, xa, qa, xb, qb, η)

    m = 1.0
    J = I(3)
    bodya = Body(m, J)
    bodyb = Body(m, J)
    bodya.state.x2[1] = xa
    bodya.state.q2[1] = qa
    bodyb.state.x2[1] = xb
    bodyb.state.q2[1] = qb
    # impulse_map_parent_jacobian_parent
    J0 = Dojo.impulse_map_jacobian(:parent, :parent, rot0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:parent, rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_parent_jacobian_child
    J0 = Dojo.impulse_map_jacobian(:parent, :child, rot0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:parent, rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_parent
    J0 = Dojo.impulse_map_jacobian(:child, :parent, rot0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:child, rot0, z[1:3], UnitQuaternion(z[4:7]..., false), xb, qb, 0)*λ,
        [xa; Dojo.vector(qa)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7

    # impulse_map_child_jacobian_child
    J0 = Dojo.impulse_map_jacobian(:child, :child, rot0, bodya, bodyb, λ)
    attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
    J1 = FiniteDiff.finite_difference_jacobian(
        z -> Dojo.impulse_map(:child, rot0, xa, qa, z[1:3], UnitQuaternion(z[4:7]..., false), 0)*λ,
        [xb; Dojo.vector(qb)]
        ) * attjac
    @test norm(J0 - J1, Inf) < 1e-7
end
