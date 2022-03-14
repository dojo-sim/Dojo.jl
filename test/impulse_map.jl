################################################################################
# Displacement Jacobian
################################################################################
@testset "Displacement Jacobian" begin
    @testset "Rotational" begin
        mech = Dojo.get_mechanism(:pendulum)
        joint0 = mech.joints[1]
        rot0 = joint0.rotational
        rot0.axis_offset = rand(QuatRotation).q

        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q

        # displacement_jacobian_configuration
        X0, Q0 = Dojo.displacement_jacobian_configuration(:parent, rot0, xa, qa, xb, qb, 
            attjac=true)
        J0 = [X0 Q0]
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.displacement(rot0, z[1:3], Quaternion(z[4:7]...), xb, qb),
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # displacement_jacobian_configuration
        X0, Q0 = Dojo.displacement_jacobian_configuration(:child, rot0, xa, qa, xb, qb, 
            attjac=true)
        J0 = [X0 Q0]
        attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.displacement(rot0, xa, qa, z[1:3], Quaternion(z[4:7]...)),
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end

    @testset "Translational" begin
        mech = Dojo.get_mechanism(:slider)
        joint0 = mech.joints[1]
        tra0 = joint0.translational

        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q

        # displacement_jacobian_configuration
        X0, Q0 = Dojo.displacement_jacobian_configuration(:parent, tra0, xa, qa, xb, qb, 
            attjac=true)
        J0 = [X0 Q0]
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.displacement(tra0, z[1:3], Quaternion(z[4:7]...), xb, qb),
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # displacement_jacobian_configuration
        X0, Q0 = Dojo.displacement_jacobian_configuration(:child, tra0, xa, qa, xb, qb, 
            attjac=true)
        J0 = [X0 Q0]
        attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.displacement(tra0, xa, qa, z[1:3], Quaternion(z[4:7]...)),
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end
end

################################################################################
# Impulse transform Jacobian
################################################################################
@testset "Impulse transform Jacobian: rotational" begin
    @testset "Rotational" begin
        mech = Dojo.get_mechanism(:pendulum)
        joint0 = mech.joints[1]
        rot0 = joint0.rotational
        rot0.axis_offset = rand(QuatRotation).q

        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q

        p0 = rand(3)
        J0 = Dojo.impulse_transform_jacobian(:parent, :parent, rot0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:parent, rot0, z[1:3], Quaternion(z[4:7]...), xb, qb) * p0,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:parent, :child, rot0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:parent, rot0, xa, qa, z[1:3], Quaternion(z[4:7]...)) * p0,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:child, :parent, rot0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:child, rot0, z[1:3], Quaternion(z[4:7]...), xb, qb) * p0,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:child, :child, rot0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:child, rot0, xa, qa, z[1:3], Quaternion(z[4:7]...)) * p0,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end

    @testset "Translational" begin
        mech = Dojo.get_mechanism(:slider)
        joint0 = mech.joints[1]
        tra0 = joint0.translational

        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q

        p0 = rand(3)
        J0 = Dojo.impulse_transform_jacobian(:parent, :parent, tra0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:parent, tra0, z[1:3], Quaternion(z[4:7]...), xb, qb) * p0,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:parent, :child, tra0, xa, qa, xb, qb, p0)
        attjac = cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:parent, tra0, xa, qa, z[1:3], Quaternion(z[4:7]...)) * p0,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:child, :parent, tra0, xa, qa, xb, qb, p0)
        attjac = cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:child, tra0, z[1:3], Quaternion(z[4:7]...), xb, qb) * p0,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        J0 = Dojo.impulse_transform_jacobian(:child, :child, tra0, xa, qa, xb, qb, p0)
        attjac = cat(I(3),Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_transform(:child, tra0, xa, qa, z[1:3], Quaternion(z[4:7]...)) * p0,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end
end


################################################################################
# Impulse map Jacobian
################################################################################
@testset "Impulse map" begin
    @testset "impulse map: rotational" begin
        mech = Dojo.get_mechanism(:pendulum)
        joint0 = mech.joints[1]
        rot0 = joint0.rotational
        rot0.axis_offset = rand(QuatRotation).q


        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q
        λ = rand(2)
        η = rand(0)
        Dojo.impulse_map(:parent, rot0, xa, qa, xb, qb, η)
        Dojo.impulse_map(:child, rot0, xa, qa, xb, qb, η)

        m = 1.0
        J = I(3)
        pbody = Body(m, J)
        cbody = Body(m, J)
        pbody.state.x2 = xa
        pbody.state.q2 = qa
        cbody.state.x2 = xb
        cbody.state.q2 = qb
        # impulse_map_parent_jacobian_parent
        J0 = Dojo.impulse_map_jacobian(:parent, :parent, rot0, pbody, cbody, λ)
        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:parent, rot0, z[1:3], Quaternion(z[4:7]...), xb, qb, 0) * λ,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_parent_jacobian_child
        J0 = Dojo.impulse_map_jacobian(:parent, :child, rot0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:parent, rot0, xa, qa, z[1:3], Quaternion(z[4:7]...), 0) * λ,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_child_jacobian_parent
        J0 = Dojo.impulse_map_jacobian(:child, :parent, rot0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:child, rot0, z[1:3], Quaternion(z[4:7]...), xb, qb, 0) * λ,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_child_jacobian_child
        J0 = Dojo.impulse_map_jacobian(:child, :child, rot0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:child, rot0, xa, qa, z[1:3], Quaternion(z[4:7]...), 0) * λ,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end

    @testset "Translational" begin
        mech = get_mechanism(:pendulum)
        joint0 = mech.joints[1]
        tra0 = joint0.translational

        xa = rand(3)
        qa = rand(QuatRotation).q
        xb = rand(3)
        qb = rand(QuatRotation).q
        λ = rand(3)
        η = rand(0)
        Dojo.impulse_map(:parent, tra0, xa, qa, xb, qb, η)
        Dojo.impulse_map(:child, tra0, xa, qa, xb, qb, η)

        m = 1.0
        J = I(3)
        pbody = Body(m, J)
        cbody = Body(m, J)
        pbody.state.x2 = xa
        pbody.state.q2 = qa
        cbody.state.x2 = xb
        cbody.state.q2 = qb
        # impulse_map_parent_jacobian_parent
        J0 = Dojo.impulse_map_jacobian(:parent, :parent, tra0, pbody, cbody, λ)

        attjac = cat(I(3),Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:parent, tra0, z[1:3], Quaternion(z[4:7]...), xb, qb, 0) * λ,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_parent_jacobian_child
        J0 = Dojo.impulse_map_jacobian(:parent, :child, tra0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:parent, tra0, xa, qa, z[1:3], Quaternion(z[4:7]...), 0) * λ,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_child_jacobian_parent
        J0 = Dojo.impulse_map_jacobian(:child, :parent, tra0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qa), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:child, tra0, z[1:3], Quaternion(z[4:7]...), xb, qb, 0) * λ,
            [xa; Dojo.vector(qa)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7

        # impulse_map_child_jacobian_child
        J0 = Dojo.impulse_map_jacobian(:child, :child, tra0, pbody, cbody, λ)
        attjac = cat(I(3), Dojo.LVᵀmat(qb), dims=(1,2))
        J1 = ForwardDiff.jacobian(
            z -> Dojo.impulse_map(:child, tra0, xa, qa, z[1:3], Quaternion(z[4:7]...), 0) * λ,
            [xb; Dojo.vector(qb)]
            ) * attjac
        @test norm(J0 - J1, Inf) < 1.0e-7
    end
end
