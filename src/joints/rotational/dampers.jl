###############################################################################
# Damper Force
###############################################################################

function damper_force(joint::Rotational{T}, qa::Quaternion, ϕa::AbstractVector,
        qb::Quaternion, ϕb::AbstractVector, 
        timestep; 
        unitary::Bool=false) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    axis_offset = joint.axis_offset
    force = 2 * Aᵀ * A * (vector_rotate(ϕb, inv(axis_offset) * inv(qa) * qb) - vector_rotate(ϕa, inv(axis_offset))) # in offset frame
    unitary && (force *= joint.damper) # Currently assumes same damper constant in all directions
    return force # in the offset frame
end

function damper_force(relative::Symbol, joint::Rotational{T}, qa::Quaternion, ϕa::AbstractVector,
        qb::Quaternion, ϕb::AbstractVector, timestep; 
        rotate::Bool=true, 
        unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    axis_offset = joint.axis_offset

    if relative == :parent
        velocity = A * (vector_rotate(ϕb, qa \ qb / axis_offset) - vector_rotate(ϕa, inv(axis_offset))) # in offset frame
        force = 2 * Aᵀ * A * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
        rotate && (force = vector_rotate(force, axis_offset)) # rotate back to frame a
        return [szeros(T, 3); force]
    elseif relative == :child 
        velocity = A * (vector_rotate(ϕb, qa \ qb / axis_offset) - vector_rotate(ϕa, inv(axis_offset))) # in offset frame
        force = - 2 * Aᵀ * A * damper * Aᵀ * velocity # currently assumes same damper constant in all directions
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
        timestep::T; attjac::Bool = true) where T

    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset
    X = szeros(T, 3, 3)
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]

    if relative == :parent 
        if jacobian == :parent 
            Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(qb * inv(axis_offset)) * Tmat()
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(inv(axis_offset)) * Lmat(inv(qa))
            attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child 
        if jacobian == :parent 
            Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(qb * inv(axis_offset)) * Tmat()
            Q += ∂vector_rotate∂q(force, inv(qb) * qa * axis_offset) * Rmat(axis_offset) * Lmat(inv(qb))
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂q(ϕb, qa \ qb / axis_offset) * Rmat(inv(axis_offset)) * Lmat(inv(qa))
            Q += ∂vector_rotate∂q(force, inv(qb) * qa * axis_offset) * Rmat(qa * axis_offset) * Tmat()
            attjac && (Q *= LVᵀmat(qb))
        end
    end 
    return timestep * [Z; X Q]
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, 
    joint::Rotational, pbody::Node, cbody::Node, 
    timestep::T) where T

    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, va, qa, ϕa = current_configuration_velocity(pbody.state)
    xb, vb, qb, ϕb = current_configuration_velocity(cbody.state)
    axis_offset = joint.axis_offset
    V = szeros(T, 3, 3)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]
    if relative == :parent 
        if jacobian == :parent 
            Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * -1.0 * ∂vector_rotate∂p(ϕa, inv(axis_offset))
        elseif jacobian == :child 
            Ω = ∂vector_rotate∂p(force, axis_offset) * 2 * joint.damper * Aᵀ * A * ∂vector_rotate∂p(ϕb, qa \ qb / axis_offset)
        end
    elseif relative == :child 
        if jacobian == :parent 
            Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * -1.0 * ∂vector_rotate∂p(ϕa, inv(axis_offset))
        elseif jacobian == :child 
            Ω = ∂vector_rotate∂p(force, inv(qb) * qa * axis_offset) * -2 * joint.damper * Aᵀ * A * ∂vector_rotate∂p(ϕb, qa \ qb / axis_offset)
        end
    end 
    return timestep * [szeros(T, 3, 6); V Ω]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T) where T = szeros(T, 6, 6)
