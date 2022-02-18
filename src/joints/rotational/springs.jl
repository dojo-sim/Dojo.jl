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
    spring_impulses(relative, joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., timestep, unitary=unitary)
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

function spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, 
    joint::Rotational, body1::Node, body2::Node,
    timestep::T; attjac::Bool = true) where T

    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    X = szeros(T, 3, 3)
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)

    if relative == :parent 
        if jacobian == :parent 
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:parent, joint, xa, UnitQuaternion(q..., false), xb, qb, timestep)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:parent, joint, xa, qa, xb, UnitQuaternion(q..., false), timestep)[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child 
        if jacobian == :parent 
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:child, joint, xa, UnitQuaternion(q..., false), xb, qb, timestep)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = FiniteDiff.finite_difference_jacobian(q -> spring_impulses(:child, joint, xa, qa, xb, UnitQuaternion(q..., false), timestep)[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    end   
    return [Z; X Q]
end

function spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return [szeros(T, 3, 6); V Ω]
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)


