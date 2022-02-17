################################################################################
# Spring Force
################################################################################

function spring_force(joint::Translational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δmincoord = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb) # in the a frame
    input = spring * Aᵀ * Δmincoord # in the a frame

    return input
end

@inline function spring_force(relative::Symbol, 
    joint::Translational{T}, 
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; 
    unitary::Bool=false) where T

    input = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
    input = impulse_transform(relative, joint, xa, qa, xb, qb) * input

    return input
end

# @inline function spring_relative(relative::Symbol, joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
#         xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
#     input = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
#     input = impulse_transform(relative, joint, xa, qa, xb, qb) * input
#     return input
# end

spring_impulses(relative::Symbol, joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_force(relative, joint,
    current_configuration(bodya.state)..., 
    current_configuration(bodyb.state)...; 
    unitary=unitary)
spring_impulses(relative::Symbol, joint::Translational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Spring Jacobians
################################################################################

function spring_jacobian_configuration(jacobian_relative::Symbol, 
        joint::Translational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    ∇input = spring * Aᵀ * - minimal_coordinates_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb)

    return ∇input
end

@inline function spring_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
    joint::Translational{T}, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep::T; unitary::Bool=false) where T

    input = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
    ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
        spring_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb, unitary=unitary)
    ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, input)

    return timestep * ∇xq
end

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    body1::Node, body2::Node, 
    timestep::T; attjac::Bool = true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    
    # attjac && (return spring_jacobian_configuration(relative, jacobian, joint, xa, qa, xb, qb, timestep; unitary=false))

    if relative == :parent 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(xa -> spring_force(:parent, joint, xa, qa, xb, qb), xa)
            Q = FiniteDiff.finite_difference_jacobian(qa -> spring_force(:parent, joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(xb -> spring_force(:parent, joint, xa, qa, xb, qb), xb)
            Q = FiniteDiff.finite_difference_jacobian(qb -> spring_force(:parent, joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(xa -> spring_force(:child, joint, xa, qa, xb, qb), xa)
            Q = FiniteDiff.finite_difference_jacobian(qa -> spring_force(:child, joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(xb -> spring_force(:child, joint, xa, qa, xb, qb), xb)
            Q = FiniteDiff.finite_difference_jacobian(qb -> spring_force(:child, joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    end

    return timestep * [X Q]
end

function spring_jacobian_velocity(relative::Symbol, 
        jacobian_relative::Symbol,
        joint::Translational,
        body1::Node, body2::Node, 
        timestep::T) where T
    return timestep * szeros(T, 6, 6)
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true, unitary::Bool=false) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T, unitary::Bool=false) where T = szeros(T, 6, 6)



