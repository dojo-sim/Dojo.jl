################################################################################
# Spring Impulses
################################################################################

function spring_force(relative::Symbol, joint::Rotational{T}, 
        xa::AbstractVector, qa::Quaternion,
        xb::AbstractVector, qb::Quaternion; 
        rotate::Bool=true, 
        unitary::Bool=false) where T

    spring = unitary ? 1.0 : joint.spring
    distance = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb)
    force = -spring * zerodimstaticadjoint(nullspace_mask(joint)) * distance # force in offset frame
    
    if relative == :parent
        rotate ? (output = vector_rotate(force, joint.axis_offset)) : (output = force) # rotate into a frame
    elseif relative == :child 
        rotate ? (output = vector_rotate(-force, inv(qb) * qa  * joint.axis_offset)) : (output = -force) # rotate back to b frame
    end

    return [szeros(T, 3); output]
end

function spring_impulses(relative::Symbol, joint::Rotational, pbody::Node, cbody::Node, timestep; 
    unitary::Bool=false)
    spring_impulses(relative, joint, 
        current_configuration(pbody.state)..., 
        current_configuration(cbody.state)..., 
        timestep, unitary=unitary)
end

function spring_impulses(relative::Symbol, joint::Rotational, 
    xa::AbstractVector, qa::Quaternion, 
    xb::AbstractVector, qb::Quaternion, 
    timestep; 
    unitary::Bool=false)
    timestep * spring_force(relative, joint, xa, qa, xb, qb; unitary=unitary)
end

spring_impulses(relative::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

################################################################################
# Spring Jacobians
################################################################################

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T}, 
    xa::AbstractVector, qa::Quaternion,
    xb::AbstractVector, qb::Quaternion,
    timestep::T; 
    rotate::Bool=true, 
    unitary::Bool=false, 
    attjac=true) where T

    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    force = spring_force(relative, joint, xa, qa, xb, qb, rotate=false)
    spring = unitary ? 1.0 : joint.spring

    if relative == :parent
        J = timestep * ∂vector_rotate∂p(force[SVector{3,Int}(4,5,6)], joint.axis_offset) * spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)
    elseif relative == :child 
        X = szeros(T, 3, 3)
       
        if rotate 
            J1 = timestep * ∂vector_rotate∂p(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.axis_offset) * -spring * zerodimstaticadjoint(nullspace_mask(joint)) * minimal_coordinates_jacobian_configuration(jacobian, joint, xa, qa, xb, qb, attjac=attjac)

            if jacobian == :parent 
                Q2 = timestep * ∂vector_rotate∂q(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.axis_offset) * Rmat(joint.axis_offset) * Lmat(inv(qb))
                attjac && (Q2 *= LVᵀmat(qa))
                J2 = [X Q2]
            elseif jacobian == :child 
                Q2 = timestep * ∂vector_rotate∂q(force[SVector{3,Int}(4,5,6)], inv(qb) * qa * joint.axis_offset) * Rmat(qa * joint.axis_offset) * Tmat()
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
    joint::Rotational, pbody::Node, cbody::Node,
    timestep::T; 
    attjac::Bool = true) where T

    spring_jacobian_configuration(relative, jacobian, joint, 
        current_configuration(pbody.state)...,
        current_configuration(cbody.state)...,
        timestep; attjac=attjac)
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, pbody::Node, cbody::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational, pbody::Node, cbody::Node, timestep::T) where T = szeros(T, 6, 6)


