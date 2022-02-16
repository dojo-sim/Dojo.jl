################################################################################
# Spring Force
################################################################################

@inline function spring_force(joint::Rotational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        unitary::Bool=false) where T
    spring = unitary ? 1.0 : joint.spring
    q = spring_extension(joint, xa, qa, xb, qb)
    distance = spring_distance(joint, q)
    input = -spring * distance # force in offset frame
    input = vrotate(input, joint.qoffset) # rotate back to the a frame
    return input
end

@inline function spring_force(relative::Symbol, joint::Rotational{T}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; 
        rotate::Bool=true, unitary::Bool=false) where T

    if relative == :parent
        force = -spring_force(joint, xa, qa, xb, qb; unitary=unitary)
    elseif relative == :child 
        force = spring_force(joint, xa, qa, xb, qb; unitary=unitary)
        rotate && (force = vrotate(force, inv(qb) * qa)) # rotate back to b frame
    end
    return [szeros(T, 3); force]
end


# @inline function spring_relative(relative::Symbol, joint::Rotational{T}, xa::AbstractVector,
#         qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
#     input = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
#     input = impulse_transform(relative, joint, xa, qa, xb, qb) * input
#     return input
# end

"""
    Sinusoidal spring model:
        This spring model has no singularity.
        torque = spring * sin(θ) * r ≡≡ spring * 2w * sqrt(1 - w^2) * r == spring * 2w * [x,y,z]
        q = [cos(θ/2); r sin(θ/2)] = [w, x, y, z]
    Linear spring model:
        This spring model has a singularity at θ = π and θ = -π. The spring torque changes direction at these points.
        We need to becareful to avoid non-differentiable point especially for θ = 0. We use sinc!
        torque = spring * θ * r == spring * θ/sin(θ/2) * sin(θ/2)r == spring * θ/sin(θ/2) * [x,y,z] == spring * 1/sinc(θ/2π) * [x,y,z] == spring * 1/sinc(atan(||r||_2, w)/π) * [x,y,z]
        q = [cos(θ/2); r sin(θ/2)] = [w, x, y, z]
"""
@inline function spring_distance(joint::Rotational, q::UnitQuaternion)
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    if joint.spring_type == :sinusoidal
        return Aᵀ*A * 2 * q.w * Vmat(q)
    elseif joint.spring_type == :linear
        nr = norm(Vmat(q))
        return Aᵀ*A * 1/sinc(atan(nr, q.w)/π) * Vmat(q)
    else
        error("unknown spring model")
    end
end

"""
    Sinusoidal spring model:
        This spring model has no singularity.
        torque = spring * sin(θ) * r ≡≡ spring * 2w * sqrt(1 - w^2) * r == spring * 2w * [x,y,z]
        energy = spring * [1 - cos(θ)] == spring * 2 * (1 - cos(θ/2)^2) == spring * 2(1 - w^2)
        q = [cos(θ/2); r sin(θ/2)] = [w, x, y, z]
    Linear spring model:
        This spring model has a singularity at θ = π and θ = -π. The spring torque changes direction at these points.
        torque = spring * θ * r == spring * θ/sin(θ/2) * sin(θ/2)r == spring * θ/sin(θ/2) * [x,y,z] == spring * 1/sinc(θ/2π) * [x,y,z] == spring * 1/sinc(atan(||r||_2, w)/π) * [x,y,z]
        energy = 1/2 * spring * θ^2 = 1/2 * spring * (2 * atan(||r||_2, w))^2
        q = [cos(θ/2); r sin(θ/2)] = [w, x, y, z]
"""
@inline function energy(joint::Rotational, q::UnitQuaternion)
    if joint.spring_type == :sinusoidal
        return joint.spring * 2 * (1 - q.w^2)
    elseif joint.spring_type == :linear
        θ = AngleAxis(q).theta
        return 0.5 * joint.spring * θ^2
    else
        error("unknown spring model")
    end
end

@inline function spring_qoffset(joint::Rotational)
    # We need to convert the joint.spring_offset which is some minimal coordinate
    # representation of the joint (i.e. axis angle). For us it's the axis-angle representation.
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    aa = Aᵀ * joint.spring_offset # axis angle
    spring_qoffset = axis_angle_to_quaternion(aa)
    return spring_qoffset
end

function spring_extension(joint::Rotational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    q = inv(spring_qoffset(joint)) * displacement(joint, xa, qa, xb, qb, vmat=false)
    return q
end

spring_force(relative::Symbol, joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_force(relative, joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)...; unitary=unitary)
spring_force(relative::Symbol, joint::Rotational{T,3}, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) where T = szeros(T, 6)

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
            Q = FiniteDiff.finite_difference_jacobian(qa -> spring_force(:parent, joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = FiniteDiff.finite_difference_jacobian(qb -> spring_force(:parent, joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    elseif relative == :child 
        if jacobian == :parent 
            Q = FiniteDiff.finite_difference_jacobian(qa -> spring_force(:child, joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
            attjac && (Q *= LVᵀmat(qa))
        elseif jacobian == :child 
            Q = FiniteDiff.finite_difference_jacobian(qb -> spring_force(:child, joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
            attjac && (Q *= LVᵀmat(qb))
        end
    end   
    return timestep * [Z; X Q]
end

function spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return timestep * [szeros(T, 3, 6); V Ω]
end

spring_jacobian_configuration(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_jacobian_velocity(relative::Symbol, jacobian::Symbol, joint::Rotational{T,3}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)


