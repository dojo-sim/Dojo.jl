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
    @show force
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


# ################################################################################
# # Damper Force
# ################################################################################
#
# function damper_force(joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#         timestep;
#         unitary::Bool=false) where T
#
#     damper = unitary ? 1.0 : joint.damper
#     Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
#     input = damper * Aᵀ * -minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
#     return input
# end
#
# function damper_force(relative::Symbol, joint::Rotational{T},
#     xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#     xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#     timestep;
#     unitary::Bool=false) where T
#
#     input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary) # in the a frame
#     inputa = impulse_transform(relative, joint, xa, qa, xb, qb) * input
#     # input = SVector{3}(0,0,-0.1)
#     # input = SVector{3}(1.0,0,0)
#     # Q = I(3)
#     # X = zeros(3,3)
#     # inputa = [X Q]' * input
#     # inputa = [szeros(T,3); vector_rotate(vector_rotate(τ, qa),inv(qb))]
#     # @show input
#     # @show inputa
#     # println("imp ", scn.(impulse_transform(relative, joint, xa, qa, xb, qb)[4:6,:]))
#     # @show xa, xb
#     # X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=true)
#     # println("Q ", scn.(Q))
#     return inputa
# end
#
# damper_impulses(relative::Symbol, joint::Rotational, pbody::Node, cbody::Node, timestep;
#     unitary::Bool=false) =
#     timestep * damper_force(relative, joint,
#         current_configuration_velocity(pbody.state)...,
#         current_configuration_velocity(cbody.state)...,
#         timestep;
#         unitary=unitary)
# damper_impulses(relative::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)
















################################################################################
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




# ################################################################################
# # Damper Jacobian
# ################################################################################
#
# function damper_force_jacobian_configuration(jacobian_relative::Symbol,
#         joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#         timestep;
#         unitary::Bool=false) where T
#
#     damper = unitary ? 1.0 : joint.damper
#     Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
#     ∇input = damper * Aᵀ * -minimal_velocities_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
#
#     return ∇input
# end
#
# function damper_force_jacobian_velocity(jacobian_relative::Symbol,
#     joint::Rotational{T},
#     xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#     xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#     timestep;
#     unitary::Bool=false) where T
#
#     damper = unitary ? 1.0 : joint.damper
#     Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
#     ∇input = damper * Aᵀ * -minimal_velocities_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
#
#     return ∇input
# end
#
# function damper_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
#         joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#         timestep::T;
#         unitary::Bool=false) where T
#
#     input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
#
#     ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
#         damper_force_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
#     ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, input)
#
#     return timestep * ∇xq
# end
#
# function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol,
#     joint::Rotational,
#     pbody::Node, cbody::Node,
#     timestep::T;
#     attjac::Bool = true) where T
#
#     return damper_jacobian_configuration(relative, jacobian,
#         joint,
#         current_configuration_velocity(pbody.state)...,
#         current_configuration_velocity(cbody.state)...,
#         timestep; unitary=false)
# end
#
# function damper_jacobian_velocity(relative::Symbol, jacobian_relative::Symbol,
#         joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::Quaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::Quaternion, ωb::AbstractVector,
#         timestep::T;
#         unitary::Bool=false) where T
#
#     ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
#         damper_force_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
#
#     return timestep * ∇xq
# end
#
# function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol,
#     joint::Rotational,
#     pbody::Node, cbody::Node,
#     timestep::T, attjac::Bool=true) where T
#
#     return damper_jacobian_velocity(relative, jacobian,
#         joint,
#         current_configuration_velocity(pbody.state)...,
#         current_configuration_velocity(cbody.state)...,
#         timestep; unitary=false)
# end
#
# damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true, unitary::Bool=false) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
# damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T, unitary::Bool=false) where T = szeros(T, 6, 6)
