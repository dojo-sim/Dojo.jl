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

@inline function damper_force(relative::Symbol, joint::Translational{T}, 
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
    xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
    timestep; unitary::Bool=false) where T

    input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary) # in the a frame
    inputa = impulse_transform(relative, joint, xa, qa, xb, qb) * input

    return inputa
end

damper_impulses(relative::Symbol, joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_force(relative, joint, 
        current_configuration_velocity(bodya.state)...,
        current_configuration_velocity(bodyb.state)..., 
        timestep; unitary=unitary)
damper_impulses(relative::Symbol, joint::Translational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Damper Jacobian
################################################################################

function damper_jacobian_configuration(jacobian_relative::Symbol, 
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
        timestep; unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇input = damper * Aᵀ * -minimal_velocities_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
    
    return ∇input
end

function damper_jacobian_velocity(jacobian_relative::Symbol, 
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector, 
        timestep; unitary::Bool=false) where T

    damper = unitary ? 1.0 : joint.damper
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇input = damper * Aᵀ * -minimal_velocities_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep) # in the a frame
    
    return ∇input
end

@inline function damper_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector, 
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T

    input = damper_force(joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_jacobian_configuration(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)
    ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, input)

    return timestep * ∇xq
end

function damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    body1::Node, body2::Node, 
    timestep::T; attjac::Bool = true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]

    # attjac && (return damper_jacobian_configuration(relative, jacobian, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    if relative == :parent 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(xa -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), xa)
            Q = FiniteDiff.finite_difference_jacobian(qa -> damper_force(:parent, joint, xa, va, UnitQuaternion(qa..., false), ωa, xb, vb, qb, ωb, timestep), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
            return timestep * [X Q]
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(xb -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), xb)
            Q = FiniteDiff.finite_difference_jacobian(qb -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, UnitQuaternion(qb..., false), ωb, timestep), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
            return timestep * [X Q]
        end
    elseif relative == :child 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(xa -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), xa)
            Q = FiniteDiff.finite_difference_jacobian(qa -> damper_force(:child, joint, xa, va, UnitQuaternion(qa..., false), ωa, xb, vb, qb, ωb, timestep), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
            return timestep * [X Q]
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(xb -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), xb)
            Q = FiniteDiff.finite_difference_jacobian(qb -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, UnitQuaternion(qb..., false), ωb, timestep), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
            return timestep * [X Q]
        end
    end
end

@inline function damper_jacobian_velocity(relative::Symbol, jacobian_relative::Symbol,
        joint::Translational{T}, 
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector, 
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector,
        timestep::T; unitary::Bool=false) where T

    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        damper_jacobian_velocity(jacobian_relative, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep, unitary=unitary)

    return timestep * ∇xq
end

function damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    body1::Node, body2::Node, 
    timestep::T, attjac::Bool=true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]

    # attjac && (return damper_jacobian_velocity(relative, jacobian, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep; unitary=false))

    if relative == :parent 
        if jacobian == :parent 
            V = FiniteDiff.finite_difference_jacobian(va -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), va)
            Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), ωa)
        elseif jacobian == :child 
            V = FiniteDiff.finite_difference_jacobian(vb -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), vb)
            Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_force(:parent, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), ωb)
        end
    elseif relative == :child 
        if jacobian == :parent 
            V = FiniteDiff.finite_difference_jacobian(va -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), va)
            Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), ωa)
        elseif jacobian == :child 
            V = FiniteDiff.finite_difference_jacobian(vb -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), vb)
            Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_force(:child, joint, xa, va, qa, ωa, xb, vb, qb, ωb, timestep), ωb)
        end
    end
    return timestep * [V Ω]
end

damper_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true, unitary::Bool=false) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T, unitary::Bool=false) where T = szeros(T, 6, 6)
