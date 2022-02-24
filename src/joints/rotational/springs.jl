################################################################################
# Spring Impulses
################################################################################

function spring_force(relative::Symbol, joint::Rotational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        rotate::Bool=true, unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    distance = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb)
    force = -spring * zerodimstaticadjoint(nullspace_mask(joint)) * distance # force in offset frame
    
    if relative == :parent
        rotate ? (output = vrotate(force, joint.qoffset)) : (output = force) # rotate into a frame
    elseif relative == :child 
        rotate ? (output = vrotate(-force, inv(qb) * qa  * joint.qoffset)) : (output = -force) # rotate back to b frame
    end

    return [szeros(T, 3); output]
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

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T}, 
    xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion,
    timestep::T; 
    rotate::Bool=true, unitary::Bool=false, attjac=true) where T

    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    force = spring_force(relative, joint, xa, qa, xb, qb, rotate=false)
    spring = unitary ? 1.0 : joint.spring

    if relative == :parent
        J = timestep * ∂vrotate∂p(force[SVector{3,Int}(4,5,6)], joint.qoffset) * spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)
    elseif relative == :child 
        X = szeros(T, 3, 3)
       
        if rotate 
            J1 = timestep * ∂vrotate∂p(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.qoffset) * -spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)

            if jacobian == :parent 
                Q2 = timestep * ∂vrotate∂q(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.qoffset) * Rmat(joint.qoffset) * Lmat(inv(qb))
                attjac && (Q2 *= LVᵀmat(qa))
                J2 = [X Q2]
            elseif jacobian == :child 
                Q2 = timestep * ∂vrotate∂q(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.qoffset) * Rmat(qa * joint.qoffset) * Tmat()
                attjac && (Q2 *= LVᵀmat(qb))
                J2 = [X Q2]
            end
            J = J1 + J2
        else 
            J = timestep * -spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)
        end
    end

    return [Z; J]
end

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Rotational, body1::Node, body2::Node,
    timestep::T; attjac::Bool = true) where T

    spring_jacobian_configuration(relative, jacobian, joint, 
        current_configuration(body1.state)...,
        current_configuration(body2.state)...,
        timestep; attjac=attjac)
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)


