@testset "rotational damper jacobian" begin
    joint_types = [
        :Fixed,
        :Prismatic,
        :Planar,
        :FixedOrientation,
        :Revolute,
        :Cylindrical,
        :PlanarAxis,
        :FreeRevolute,
        :Orbital,
        :PrismaticOrbital,
        :PlanarOrbital,
        :FreeOrbital,
        :Spherical,
        :CylindricalFree,
        :PlanarFree
        ]
    for joint_type in joint_types
        mech = DojoEnvironments.get_snake(gravity=0.0, num_bodies=2, damper=0.3, joint_type=joint_type)
        DojoEnvironments.initialize_snake!(mech)
        function ctrl!(m,k)
            Dojo.set_input!(m, 0.01*m.timestep*ones(Dojo.minimal_dimension(m)))
        end
        storage = Dojo.simulate!(mech, 1.0, ctrl!)
        # Dojo.visualize(mech, storage, vis=vis)


        rot0 = mech.joints[2].rotational
        timestep0 = mech.timestep
        pbody0 = mech.bodies[1]
        cbody0 = mech.bodies[2]
        xa0, va0, qa0, ϕa0 = Dojo.current_configuration_velocity(pbody0.state)
        xb0, vb0, qb0, ϕb0 = Dojo.current_configuration_velocity(cbody0.state)

        # Configuration
        J0 = Dojo.damper_jacobian_configuration(:parent, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            xq -> timestep0 * Dojo.damper_force(:parent, rot0,
                Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_configuration(:parent, :child, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            xq -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0,
                Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
            [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_configuration(:child, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            xq -> timestep0 * Dojo.damper_force(:child, rot0,
                Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_configuration(:child, :child, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            xq -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0,
                Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
            [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        # Velocity
        J0 = Dojo.damper_jacobian_velocity(:parent, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0,
                vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [va0; ϕa0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_velocity(:parent, :child, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, qb0,
                vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
            [vb0; ϕb0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_velocity(:child, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0,
                vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [va0; ϕa0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8

        J0 = Dojo.damper_jacobian_velocity(:child, :child, rot0, pbody0, cbody0, timestep0)
        J1 = ForwardDiff.jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, qb0,
                vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
            [vb0; ϕb0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1.0e-8
    end
end
