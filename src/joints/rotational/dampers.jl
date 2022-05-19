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

    velocity = minimal_velocities(joint, szeros(3), szeros(3), qa, ϕa, szeros(3), szeros(3), qb, ϕb, timestep)
    if relative == :parent
        force = 1.0000 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, axis_offset)) # rotate back to frame a
        return [szeros(T, 3); force]
    elseif relative == :child
        force = - 1.0000 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
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
        ) where T

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset
    X = szeros(T, 3, 3)
    Z = szeros(T, 3, 6)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SUnitRange(4,6)]
    ∂vel = minimal_velocities_jacobian_configuration(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)[:, SUnitRange(4,6)]

    if relative == :parent
        Q = rotation_matrix(axis_offset) * 1.0000 * joint.damper * Aᵀ * ∂vel
    elseif relative == :child
        Q = rotation_matrix(inv(qb) * qa * axis_offset) * -1.0000 * joint.damper * Aᵀ * ∂vel
        if jacobian == :parent
            # Q += rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, vector_rotate(force, axis_offset), attjac=true)
            Q += rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, vector_rotate(force, axis_offset)) * LVᵀmat(qa) # ATTJAC
        elseif jacobian == :child
            # Q += ∂rotation_matrix_inv∂q(qb, vector_rotate(force, qa * axis_offset), attjac=true)
            Q += ∂rotation_matrix_inv∂q(qb, vector_rotate(force, qa * axis_offset)) * LVᵀmat(qb) # ATTJAC
        end
    end
    return timestep * [Z; X Q]
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol,
    joint::Rotational, pbody::Node, cbody::Node,
    timestep::T) where T

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SUnitRange(4,6)]
    ∂vel = minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)

    if relative == :parent
        VΩ = rotation_matrix(axis_offset) * 1.0000 * joint.damper * Aᵀ * ∂vel
    elseif relative == :child
        VΩ = rotation_matrix(inv(qb) * qa * axis_offset) * -1.0000 * joint.damper * Aᵀ * ∂vel
    end
    return timestep * [szeros(T, 3, 6); VΩ]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T) where T = szeros(T, 6, 6)



#
#
#
# using Test
# using Plots
# using FiniteDiff
# vis = Visualizer()
# open(vis)
#
#
#
# mech = Dojo.get_snake(gravity=0.00, num_bodies=2, damper=0.3, joint_type=:Revolute)
# Dojo.initialize_snake!(mech)
# function ctrl!(m,k)
#     set_input!(m, 0.01*m.timestep*ones(minimal_dimension(m)))
# end
# storage = Dojo.simulate!(mech, 1.0, ctrl!)
# Dojo.visualize(mech, storage, vis=vis)
#
#
# rot0 = mech.joints[2].rotational
# timestep0 = mech.timestep
# pbody0 = mech.bodies[1]
# cbody0 = mech.bodies[2]
# xa0, va0, qa0, ϕa0 = Dojo.current_configuration_velocity(pbody0.state)
# xb0, vb0, qb0, ϕb0 = Dojo.current_configuration_velocity(cbody0.state)
#
# # Configuration
# J0 = Dojo.damper_jacobian_configuration(:parent, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:parent, rot0,
#         Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_configuration(:parent, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0,
#         Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
#     [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_configuration(:child, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:child, rot0,
#         Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
#
# J0 = Dojo.damper_jacobian_configuration(:child, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0,
#         Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
#     [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# mech = Dojo.get_snake(gravity=0.00, num_bodies=2, damper=0.3, joint_type=:Revolute)
# Dojo.initialize_snake!(mech)
# function ctrl!(m,k)
#     set_input!(m, 0.01*m.timestep*ones(minimal_dimension(m)))
# end
# storage = Dojo.simulate!(mech, 1.0, ctrl!)
# Dojo.visualize(mech, storage, vis=vis)
#
# mech.joints[2]
# rot0 = mech.joints[2].rotational
# timestep0 = mech.timestep
# pbody0 = mech.bodies[1]
# cbody0 = mech.bodies[2]
# xa0, va0, qa0, ϕa0 = Dojo.current_configuration_velocity(pbody0.state)
# xb0, vb0, qb0, ϕb0 = Dojo.current_configuration_velocity(cbody0.state)
#
# # Configuration
# J0 = Dojo.damper_jacobian_configuration(:parent, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:parent, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_configuration(:parent, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
#     [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
# norm(J0 - J1, Inf)
# J0 + J1
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_configuration(:child, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:child, rot0, Dojo.Quaternion(xq[4:7]...,true), ϕa0, qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
# J0
# J1
#
# J0 = Dojo.damper_jacobian_configuration(:child, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, Dojo.Quaternion(xq[4:7]...,true), ϕb0, timestep0; rotate=true, unitary=false),
#     [xb0; Dojo.vector(qb0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qb0), dims=(1,2))
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
# J0
# J1
#
#
#
# timestep0 = 0.01
# xa0 = szeros(3)
# qa0 = Quaternion(1,0,0,0.0,true)
# va0 = szeros(3)
# ϕa0 = szeros(3)
# xb0 = szeros(3)
# qb0 = Quaternion(1,0,0,0.0,true)
# vb0 = szeros(3)
# ϕb0 = szeros(3)
# J0 = minimal_velocities_jacobian_velocity(:parent,
#     rot0, xa0, va0, qa0, ϕa0, xb0, vb0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     xq -> Dojo.minimal_velocities(rot0, xq[SUnitRange(1,3)], va0, Dojo.Quaternion(xq[4:7]...,true), ϕa0,
#         xb0, vb0, qb0, ϕb0, timestep0),
#     [xa0; Dojo.vector(qa0)]) * Dojo.cat(I(3), Dojo.LVᵀmat(qa0), dims=(1,2))
# norm(J0 - J1, Inf)
#
# using Rotations
#
#
#
#
#
#
#
#
#
#
#
#
# q0 = rand(QuatRotation).q
# J0 = drotation_vectordq(q0)
# J1 = FiniteDiff.finite_difference_jacobian(q -> rotation_vector(q), vector(q0))
# norm(J0 - J1, Inf)
#
# timestep0 = 0.01
# q0 = rand(QuatRotation).q
# ϕ0 = srand(3)
# J0 = rotational_integrator_jacobian_velocity(q0, -ϕ0, timestep0)
# J1 = -FiniteDiff.finite_difference_jacobian(ϕ -> vector(next_orientation(q0, -ϕ, timestep0)), ϕ0)
# norm(J0 - J1, Inf)
#
#
#
#
#
#
# function minvel(joint::Rotational, qb::Quaternion, ϕb::AbstractVector, timestep)
# 	A = nullspace_mask(joint)
# 	qb1 = next_orientation(qb, -ϕb, timestep)
#     return A * rotation_vector(qb1)
# end
#
# function minvel_jacobian_velocity(joint::Rotational{T}, qb::Quaternion,
#         ϕb::AbstractVector, timestep) where T
#
# 	A = nullspace_mask(joint)
#     qb1 = next_orientation(qb, -ϕb, timestep)
#     Ω = A * drotation_vectordq(qb1) * -rotational_integrator_jacobian_velocity(qb, -ϕb, timestep)
# 	return Ω
# end
#
# LVᵀmat(qb0) * VLᵀmat(qb0)
#
# timestep0 = 1.0
# xa0 = srand(3)
# # qa0 = rand(QuatRotation).q
# qa0 = Quaternion(1,0,0,0.0,true)
# va0 = srand(3)
# ϕa0 = srand(3)
# xb0 = srand(3)
# qb0 = rand(QuatRotation).q
# # qb0 = Quaternion(1,0,0,0.0,true)
# vb0 = srand(3)
# ϕb0 = srand(3)
# rot0
# J0 = minvel_jacobian_velocity(rot0, qb0, ϕb0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     ϕ -> Dojo.minvel(rot0, qb0, ϕ, timestep0),
#     ϕb0)
# norm(J0 - J1, Inf)
#
#
#
#
#
#
#
#
#
#
# # Velocity
# J0 = Dojo.damper_jacobian_velocity(:parent, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [va0; ϕa0])
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_velocity(:parent, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> timestep0 * Dojo.damper_force(:parent, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
#     [vb0; ϕb0])
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_velocity(:child, :parent, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, vϕ[Dojo.SUnitRange(4,6)], qb0, ϕb0, timestep0; rotate=true, unitary=false),
#     [va0; ϕa0])
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
# J0 = Dojo.damper_jacobian_velocity(:child, :child, rot0, pbody0, cbody0, timestep0)
# J1 = FiniteDiff.finite_difference_jacobian(
#     vϕ -> timestep0 * Dojo.damper_force(:child, rot0, qa0, ϕa0, qb0, vϕ[Dojo.SUnitRange(4,6)], timestep0; rotate=true, unitary=false),
#     [vb0; ϕb0])
# norm(J0 - J1, Inf)
# @test norm(J0 - J1, Inf) < 1e-6
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
