################################################################################
# Damper Force
################################################################################

function damper_force(joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
        timestep; unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    input = damper * Aᵀ * -minimal_velocities(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame

    return input
end

function damper_force(relative::Symbol, joint::Translational{T}, 
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
    xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
    timestep; unitary::Bool=false) where T

    input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary) # in the a frame
    inputa = impulse_transform(relative, joint, xa, qa, xb, qb) * input

    return inputa
end

damper_impulses(relative::Symbol, joint::Translational, pbody::Node, cbody::Node, timestep; unitary::Bool=false) =
    timestep * damper_force(relative, joint, 
        current_configuration_velocity(pbody.state)...,
        current_configuration_velocity(cbody.state)..., 
        timestep; unitary=unitary)
damper_impulses(relative::Symbol, joint::Translational{T,3}, pbody::Node, cbody::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Damper Jacobian
################################################################################

function damper_force_jacobian_configuration(jacobian_relative::Symbol, 
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
        timestep; unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇input = damper * Aᵀ * -minimal_velocities_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
    
    return ∇input
end

function damper_force_jacobian_velocity(jacobian_relative::Symbol, 
    joint::Translational{T}, 
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
    xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
    timestep; unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇input = damper * Aᵀ * -minimal_velocities_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame

    return ∇input
end

function damper_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector, 
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T

    input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
    
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_force_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
    ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, input)

    return timestep * ∇xq
end

function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    pbody::Node, cbody::Node, 
    timestep::T; attjac::Bool = true) where T

    return damper_jacobian_configuration(relative, jacobian, 
        joint, 
        current_configuration_velocity(pbody.state)..., 
        current_configuration_velocity(cbody.state)...,
        timestep; unitary=false)
end

function damper_jacobian_velocity(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector, 
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T

    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_force_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)

    return timestep * ∇xq
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    pbody::Node, cbody::Node, 
    timestep::T, attjac::Bool=true) where T

    return damper_jacobian_velocity(relative, jacobian, 
        joint, 
        current_configuration_velocity(pbody.state)...,
        current_configuration_velocity(cbody.state)..., 
        timestep; unitary=false)
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true, unitary::Bool=false) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, pbody::Node, cbody::Node, timestep::T, unitary::Bool=false) where T = szeros(T, 6, 6)
