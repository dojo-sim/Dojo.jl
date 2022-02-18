################################################################################
# Spring Impulses
################################################################################

function spring_force(relative::Symbol, joint::Translational{T}, 
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; 
    unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    distance = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb)
    force = spring  * distance

    return force
end

function spring_impulses(relative::Symbol, joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false)
    spring_impulses(relative, joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., timestep, unitary=unitary)
end

function spring_impulses(relative::Symbol, joint::Translational, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep; unitary::Bool=false)
    timestep * impulse_transform(relative, joint, xa, qa, xb, qb) * zerodimstaticadjoint(nullspace_mask(joint)) * spring_force(relative, joint, xa, qa, xb, qb; unitary=unitary)
end

spring_impulses(relative::Symbol, joint::Translational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Spring Jacobians
################################################################################

# function spring_jacobian_configuration(jacobian_relative::Symbol, 
#         joint::Translational{T}, 
#         xa::AbstractVector, qa::UnitQuaternion,
#         xb::AbstractVector, qb::UnitQuaternion; 
#         unitary::Bool=false) where T

#     spring = unitary ? 1.0 : joint.spring
#     Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
#     ∇input = spring * Aᵀ * - minimal_coordinates_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb)

#     return ∇input
# end

# @inline function spring_jacobian_configuration(relative::Symbol, jacobian_relative::Symbol,
#     joint::Translational{T}, 
#     xa::AbstractVector, qa::UnitQuaternion, 
#     xb::AbstractVector, qb::UnitQuaternion, 
#     timestep::T; unitary::Bool=false) where T

#     input = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
#     ∇xq = impulse_transform(relative, joint, xa, qa, xb, qb) *
#         spring_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb, unitary=unitary)
#     ∇xq += impulse_transform_jacobian(relative, jacobian_relative, joint, xa, qa, xb, qb, input)

#     return timestep * ∇xq
# end

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Translational, 
    body1::Node, body2::Node, 
    timestep::T; attjac::Bool = true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)

    # attjac && (return spring_jacobian_configuration(relative, jacobian, joint, xa, qa, xb, qb, timestep; unitary=false))

    if relative == :parent 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(x -> spring_impulses(:parent, joint, x, qa, xb, qb, timestep, unitary=unitary), xa)
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:parent, joint, xa, UnitQuaternion(q..., false), xb, qb, timestep, unitary=unitary), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(x -> spring_impulses(:parent, joint, xa, qa, x, qb, timestep, unitary=unitary), xb)
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:parent, joint, xa, qa, xb, UnitQuaternion(q..., false), timestep, unitary=unitary), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child 
        if jacobian == :parent 
            X = FiniteDiff.finite_difference_jacobian(x -> spring_impulses(:child, joint, x, qa, xb, qb, timestep, unitary=unitary), xa)
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:child, joint, xa, UnitQuaternion(q..., false), xb, qb, timestep, unitary=unitary), [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            X = FiniteDiff.finite_difference_jacobian(x -> spring_impulses(:child, joint, xa, qa, x, qb, timestep, unitary=unitary), xb)
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:child, joint, xa, qa, xb, UnitQuaternion(q..., false), timestep, unitary=unitary), [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    end

    return [X Q]
end

function spring_jacobian_velocity(relative::Symbol, 
    jacobian_relative::Symbol,
    joint::Translational,
    body1::Node, body2::Node, 
    timestep::T) where T
    return szeros(T, 6, 6)
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true, unitary::Bool=false) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Translational{T,3}, body1::Node, body2::Node, timestep::T, unitary::Bool=false) where T = szeros(T, 6, 6)



