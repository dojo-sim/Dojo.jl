###############################################################################
# Damper Force
###############################################################################

function damper_force(relative::Symbol, joint::Rotational{T}, qa::Quaternion, ϕa::AbstractVector,
        qb::Quaternion, ϕb::AbstractVector, timestep;
        rotate::Bool=true,
        unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    axis_offset = joint.axis_offset

    if relative == :parent
        # velocity = A * (vector_rotate(ϕb, qa \ qb / axis_offset) - vector_rotate(ϕa, inv(axis_offset))) # in offset frame
        velocity = minimal_velocities(joint, szeros(3), szeros(3), qa, ϕa, szeros(3), szeros(3), qb, ϕb, timestep)
        force = 2 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, axis_offset)) # rotate back to frame a
        return [szeros(T, 3); force]
    elseif relative == :child
        # velocity = A * (vector_rotate(ϕb, qa \ qb / axis_offset) - vector_rotate(ϕa, inv(axis_offset))) # in offset frame
        velocity = minimal_velocities(joint, szeros(3), szeros(3), qa, ϕa, szeros(3), szeros(3), qb, ϕb, timestep)
        force = - 2 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, inv(qb) * qa * axis_offset)) # rotate back to frame b
        return [szeros(T, 3); force]
    end
end

damper_impulses(relative::Symbol, joint::Rotational, pbody::Node, cbody::Node, timestep; unitary::Bool=false) =
    timestep * damper_force(relative, joint, current_configuration(pbody.state)[2], pbody.state.ϕsol[2],
    current_configuration(cbody.state)[2], cbody.state.ϕsol[2], timestep, unitary=unitary)

damper_impulses(relative::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

###############################################################################
# Damper Jacobians
################################################################################

function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol,
        joint::Rotational, pbody::Node, cbody::Node,
        timestep::T;
        # attjac::Bool = true
        ) where T

    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset
    X = szeros(T, 3, 3)
    # Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    Z = szeros(T, 3, 6)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]

    if relative == :parent
        if jacobian == :parent
            # Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(qb * inv(axis_offset)) * Tmat()
            Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * minimal_velocities_jacobian_configuration(relative, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:, SUnitRange(4,6)]
            # attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child
            # Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(inv(axis_offset)) * Lmat(inv(qa))
            Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * minimal_velocities_jacobian_configuration(relative, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:, SUnitRange(4,6)]
            # attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child
        if jacobian == :parent
            # Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(qb * inv(axis_offset)) * Tmat()
            Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * minimal_velocities_jacobian_configuration(relative, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:, SUnitRange(4,6)]
            Q += rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, vector_rotate(force, axis_offset), attjac=true)
            # Q += ∂vector_rotate∂q(force, inv(qb) * qa * axis_offset) * Rmat(axis_offset) * Lmat(inv(qb))
            # attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child
            # Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(inv(axis_offset)) * Lmat(inv(qa))
            Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * minimal_velocities_jacobian_configuration(relative, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:, SUnitRange(4,6)]
            # Q += ∂vector_rotate∂q(force, inv(qb) * qa * axis_offset) * Rmat(qa * axis_offset) * Tmat()
            Q += ∂rotation_matrix_inv∂q(qb, vector_rotate(force, qa * axis_offset), attjac=true)
            # attjac && (Q *= LVᵀmat(qb))
        end
    end
    return timestep * [Z; X Q]
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol,
    joint::Rotational, pbody::Node, cbody::Node,
    timestep::T) where T
    # @show "yo"
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset
    V = szeros(T, 3, 3)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]
    if relative == :parent
        if jacobian == :parent
            # Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * -1.0 * ∂vector_rotate∂p(ϕa, inv(axis_offset))
            Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:,SUnitRange(4,6)]
        elseif jacobian == :child
            # Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂p(ϕb, qa \ qb / axis_offset)
            Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:,SUnitRange(4,6)]
        end
    elseif relative == :child
        if jacobian == :parent
            # Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * -1.0 * ∂vector_rotate∂p(ϕa, inv(axis_offset))
            Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:,SUnitRange(4,6)]
        elseif jacobian == :child
            # Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂p(ϕb, qa \ qb / axis_offset)
            Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:,SUnitRange(4,6)]
        end
    end
    return timestep * [szeros(T, 3, 6); V Ω]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T) where T = szeros(T, 6, 6)


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
        mech = Dojo.get_npendulum(num_bodies=2, damper=0.3, joint_type=joint_type)
        Dojo.initialize_npendulum!(mech, base_angle=0.5)
        storage = Dojo.simulate!(mech, 1.0)
        # Dojo.visualize(mech, storage, vis=vis)

        rot0 = mech.joints[1].rotational
        timestep0 = mech.timestep
        pbody0 = mech.bodies[1]
        cbody0 = mech.bodies[2]
        xa0, va0, qa0, ϕa0 = Dojo.current_configuration_velocity(pbody0.state)
        xb0, vb0, qb0, ϕb0 = Dojo.current_configuration_velocity(cbody0.state)

        # Configuration
        J0 = Dojo.damper_jacobian_configuration(:parent, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            xq -> timestep0 * Dojo.damper_force(:parent, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_configuration(:parent, :child, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            xq -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
            [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_configuration(:child, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            xq -> timestep0 * Dojo.damper_force(:child, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_configuration(:child, :child, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            xq -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
            [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        # Velocity
        J0 = Dojo.damper_jacobian_velocity(:parent, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [va0; ϕa0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_velocity(:parent, :child, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
            [vb0; ϕb0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_velocity(:child, :parent, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
            [va0; ϕa0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6

        J0 = Dojo.damper_jacobian_velocity(:child, :child, rot0, pbody0, cbody0, timestep0)
        J1 = FiniteDiff.finite_difference_jacobian(
            vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
            [vb0; ϕb0])
        norm(J0 - J1, Inf)
        @test norm(J0 - J1, Inf) < 1e-6
    end
end





mech = Dojo.get_npendulum(num_bodies=2, damper=0.3, joint_type=:PlanarAxis)
Dojo.initialize_npendulum!(mech, base_angle=0.5)
storage = Dojo.simulate!(mech, 1.0)
Dojo.visualize(mech, storage, vis=vis)

rot0 = mech.joints[1].rotational
timestep0 = mech.timestep
pbody0 = mech.bodies[1]
cbody0 = mech.bodies[2]
xa0, va0, qa0, ϕa0 = Dojo.current_configuration_velocity(pbody0.state)
xb0, vb0, qb0, ϕb0 = Dojo.current_configuration_velocity(cbody0.state)

# Configuration
J0 = Dojo.damper_jacobian_configuration(:parent, :parent, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    xq -> timestep0 * Dojo.damper_force(:parent, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
    [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_configuration(:parent, :child, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    xq -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
    [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_configuration(:child, :parent, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    xq -> timestep0 * Dojo.damper_force(:child, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
    [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_configuration(:child, :child, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    xq -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
    [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

# Velocity
J0 = Dojo.damper_jacobian_velocity(:parent, :parent, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
    [va0; ϕa0])
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_velocity(:parent, :child, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
    [vb0; ϕb0])
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_velocity(:child, :parent, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
    [va0; ϕa0])
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6

J0 = Dojo.damper_jacobian_velocity(:child, :child, rot0, pbody0, cbody0, timestep0)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
    [vb0; ϕb0])
norm(J0 - J1, Inf)
@test norm(J0 - J1, Inf) < 1e-6
