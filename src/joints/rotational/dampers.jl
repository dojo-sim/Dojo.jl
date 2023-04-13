###############################################################################
# Damper Force
###############################################################################
function damper_force(relative::Symbol, joint::Rotational{T}, qa::Quaternion, ωa::AbstractVector,
        qb::Quaternion, ωb::AbstractVector, timestep;
        rotate::Bool=true,
        unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    orientation_offset = joint.orientation_offset

    velocity = minimal_velocities(joint, szeros(3), szeros(3), qa, ωa, szeros(3), szeros(3), qb, ωb, timestep)
    if relative == :parent
        force = 1.0 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, orientation_offset)) # rotate back to frame a
        return [szeros(T, 3); force]
    elseif relative == :child
        force = - 1.0 * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, inv(qb) * qa * orientation_offset)) # rotate back to frame b
        return [szeros(T, 3); force]
    end
end

damper_impulses(relative::Symbol, joint::Rotational, pbody::Node, cbody::Node, timestep; unitary::Bool=false) =
    timestep * damper_force(relative, joint, current_configuration(pbody.state)[2], pbody.state.ωsol[2],
    current_configuration(cbody.state)[2], cbody.state.ωsol[2], timestep, unitary=unitary)

damper_impulses(relative::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

###############################################################################
# Damper Jacobians
################################################################################

function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol,
        joint::Rotational, pbody::Node, cbody::Node,
        timestep::T;
        ) where T

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xa, va, qa, ωa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ωb = current_configuration_velocity(cbody.state)
    orientation_offset = joint.orientation_offset
    X = szeros(T, 3, 3)
    Z = szeros(T, 3, 6)

    force = damper_force(relative, joint, qa, ωa, qb, ωb, timestep; rotate = false)[SUnitRange(4,6)]
    ∂vel = minimal_velocities_jacobian_configuration(jacobian, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)[:, SUnitRange(4,6)]

    if relative == :parent
        Q = rotation_matrix(orientation_offset) * 1.0000 * joint.damper * Aᵀ * ∂vel
    elseif relative == :child
        Q = rotation_matrix(inv(qb) * qa * orientation_offset) * -1.0000 * joint.damper * Aᵀ * ∂vel
        if jacobian == :parent
            # Q += rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, vector_rotate(force, orientation_offset), attjac=true)
            Q += rotation_matrix(inv(qb)) * ∂rotation_matrix∂q(qa, vector_rotate(force, orientation_offset)) * LVᵀmat(qa) # ATTJAC
        elseif jacobian == :child
            # Q += ∂rotation_matrix_inv∂q(qb, vector_rotate(force, qa * orientation_offset), attjac=true)
            Q += ∂rotation_matrix_inv∂q(qb, vector_rotate(force, qa * orientation_offset)) * LVᵀmat(qb) # ATTJAC
        end
    end
    return timestep * [Z; X Q]
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol,
    joint::Rotational, pbody::Node, cbody::Node,
    timestep::T) where T

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xa, va, qa, ωa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ωb = current_configuration_velocity(cbody.state)
    orientation_offset = joint.orientation_offset

    force = damper_force(relative, joint, qa, ωa, qb, ωb, timestep; rotate = false)[SUnitRange(4,6)]
    ∂vel = minimal_velocities_jacobian_velocity(jacobian, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep)

    if relative == :parent
        VΩ = rotation_matrix(orientation_offset) * 1.0000 * joint.damper * Aᵀ * ∂vel
    elseif relative == :child
        VΩ = rotation_matrix(inv(qb) * qa * orientation_offset) * -1.0000 * joint.damper * Aᵀ * ∂vel
    end
    return timestep * [szeros(T, 3, 6); VΩ]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T) where T = szeros(T, 6, 6)
