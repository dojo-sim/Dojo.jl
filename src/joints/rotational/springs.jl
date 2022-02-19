################################################################################
# Spring Impulses
################################################################################

@inline function spring_force(relative::Symbol, joint::Rotational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        rotate::Bool=true, unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    distance = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb)
    force = -spring * zerodimstaticadjoint(nullspace_mask(joint)) * distance # force in offset frame
    
    if relative == :parent
        nothing
    elseif relative == :child 
        rotate && (force = vrotate(-force, inv(qb) * qa)) # rotate back to b frame
    end

    return [szeros(T, 3); force]
end

function spring_impulses(relative::Symbol, joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false)
    spring_impulses(relative, joint, 
        current_configuration(bodya.state)..., 
        current_configuration(bodyb.state)..., 
        timestep, unitary=unitary)
end

function spring_impulses(relative::Symbol, joint::Rotational, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion, 
    timestep; unitary::Bool=false)
    timestep * spring_force(relative, joint, xa, qa, xb, qb; unitary=unitary)
end

spring_impulses(relative::Symbol, joint::Rotational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Spring Jacobians
################################################################################

@inline function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T}, 
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion,
    timestep::T; 
    rotate::Bool=true, unitary::Bool=false, attjac=true) where T

    spring = unitary ? 1.0 : joint.spring
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)

    if relative == :parent
        J = spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)
    elseif relative == :child 
        X = szeros(T, 3, 3)

        input = -spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates(joint, xa, qa, xb, qb)
        input_jacobian = -spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)
       
        if rotate 
            J1 = ∂vrotate∂p(input, inv(qb) * qa) * input_jacobian
            if jacobian == :parent
                Q2 = ∂vrotate∂q(input, inv(qb) * qa) * Lmat(inv(qb))
                attjac && (Q2 *= LVᵀmat(qa))
            elseif jacobian == :child 
                Q2 = ∂vrotate∂q(input, inv(qb) * qa) * Rmat(qa) * Tmat()
                attjac && (Q2 *= LVᵀmat(qb))
            end
            J = J1 + [X Q2] 
        else
            J = input_jacobian 
        end
    end

    return timestep * [Z; J]
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)


