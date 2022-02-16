###############################################################################
# Damper Force
###############################################################################

@inline function damper_force(joint::Rotational{T}, qa::UnitQuaternion, ϕa::AbstractVector,
        qb::UnitQuaternion, ϕb::AbstractVector, timestep; unitary::Bool=false) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    force = 2Aᵀ * A * (vrotate(ϕb, inv(qoffset) * inv(qa) * qb) - vrotate(ϕa, inv(qoffset))) # in offset frame
    unitary && (force *= joint.damper) # Currently assumes same damper constant in all directions
    return force # in the offset frame
end

@inline function damper_force(relative::Symbol, joint::Rotational{T}, qa::UnitQuaternion, ϕa::AbstractVector,
        qb::UnitQuaternion, ϕb::AbstractVector, timestep; rotate::Bool=true, unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset

    if relative == :parent
        velocity = A * (vrotate(ϕb, qa \ qb / qoffset) - vrotate(ϕa, inv(qoffset))) # in offset frame
        force = 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
        rotate && (force = vrotate(force, qoffset)) # rotate back to frame a
        return [szeros(T, 3); force]
    elseif relative == :child 
        velocity = A * (vrotate(ϕb, qa \ qb / qoffset) - vrotate(ϕa, inv(qoffset))) # in offset frame
        force = - 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
        rotate && (force = vrotate(force, inv(qb) * qa * qoffset)) # rotate back to frame b
        return [szeros(T, 3); force]
    end
end

damper_force(relative::Symbol, joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_force(relative, joint, current_configuration(bodya.state)[2], bodya.state.ϕsol[2],
    current_configuration(bodyb.state)[2], bodyb.state.ϕsol[2], timestep, unitary=unitary)

damper_force(relative::Symbol, joint::Rotational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

###############################################################################
# Damper Jacobians
################################################################################

function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
        joint::Rotational, body1::Node, body2::Node,
        timestep::T; attjac::Bool = true) where T

    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ϕa = current_configuration_velocity(body1.state)
    _, _, _, ϕb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    X = szeros(T, 3, 3)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]

    if relative == :parent 
        if jacobian == :parent 
            Q = ∂vrotate∂p(force, qoffset) * 2 * joint.damper * Aᵀ * A * ∂vrotate∂q(ϕb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
            attjac && (Q *= LVᵀmat(qa))
            Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
        elseif jacobian == :child 
            Q = ∂vrotate∂p(force, qoffset) * 2 * joint.damper * Aᵀ * A * ∂vrotate∂q(ϕb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
            attjac && (Q *= LVᵀmat(qb))
            Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
        end
    elseif relative == :child 
        if jacobian == :parent 
            Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * joint.damper * Aᵀ * A * ∂vrotate∂q(ϕb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
            Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb))
            attjac && (Q *= LVᵀmat(qa))
            Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
        elseif jacobian == :child 
            Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * joint.damper * Aᵀ * A * ∂vrotate∂q(ϕb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
            Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat()
            attjac && (Q *= LVᵀmat(qb))
            Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
        end
    end 
    return timestep * [Z; X Q]
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, 
    joint::Rotational, body1::Node, body2::Node, 
    timestep::T) where T

    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ϕa = current_configuration_velocity(body1.state)
    _, _, _, ϕb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    V = szeros(T, 3, 3)

    force = damper_force(relative, joint, qa, ϕa, qb, ϕb, timestep; rotate = false)[SVector{3,Int}(4,5,6)]
    if relative == :parent 
        if jacobian == :parent 
            Ω = ∂vrotate∂p(force, qoffset) * 2 * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ϕa, inv(qoffset))
        elseif jacobian == :child 
            Ω = ∂vrotate∂p(force, qoffset) * 2 * joint.damper * Aᵀ * A * ∂vrotate∂p(ϕb, qa \ qb / qoffset)
        end
    elseif relative == :child 
        if jacobian == :parent 
            Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ϕa, inv(qoffset))
        elseif jacobian == :child 
            Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * joint.damper * Aᵀ * A * ∂vrotate∂p(ϕb, qa \ qb / qoffset)
        end
    end 
    return timestep * [szeros(T, 3, 6); V Ω]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
