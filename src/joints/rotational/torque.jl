@inline function rotation_error(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion; qoff::UnitQuaternion = spring_qoffset(joint))
    q = qa \ qb / joint.qoffset / qoff
end

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
    q = inv(spring_qoffset(joint)) * orientation_error(joint, xa, qa, xb, qb)
    return q
end

spring_parent(joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_parent(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
spring_child(joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_child(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
damper_parent(joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_parent(joint, current_configuration(bodya.state)[2], bodya.state.ϕsol[2], current_configuration(bodyb.state)[2], bodyb.state.ϕsol[2], unitary=unitary)
damper_child(joint::Rotational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_child(joint, current_configuration(bodya.state)[2], bodya.state.ϕsol[2], current_configuration(bodyb.state)[2], bodyb.state.ϕsol[2], unitary=unitary)

spring_parent(joint::Rotational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
spring_child(joint::Rotational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_parent(joint::Rotational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_child(joint::Rotational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)

spring_parent(joint::Rotational3{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)
spring_child(joint::Rotational3{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)


################################################################################
# Spring Force
################################################################################
@inline function spring_force(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    spring = unitary ? 1.0 : joint.spring
    q = spring_extension(joint, xa, qa, xb, qb)
    distance = spring_distance(joint, q)
    Fτ = -spring * distance # force in offset frame
    Fτ = vrotate(Fτ, joint.qoffset) # rotate back to the a frame
    return Fτ
end

@inline function spring_relative(relative::Symbol, joint::Rotational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; unitary::Bool=false) where T
    Fτ = spring_force(joint, xa, qa, xb, qb, unitary=unitary)
    Fτ = impulse_transform(relative, joint, xa, qa, xb, qb) * Fτ
    return Fτ
end

@inline function spring_parent(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    # spring = unitary ? 1.0 : joint.spring
    # # q = rotation_error(joint, qa, qb, qoff = spring_qoffset(joint))
    # q = spring_extension(joint, xa, qa, xb, qb)
    # distance = spring_distance(joint, q)
    # force = spring * distance # force in offset frame
    # rotate && (force = vrotate(force, joint.qoffset)) # rotate back to a frame
    force = -spring_force(joint, xa, qa, xb, qb; unitary=unitary)
    # force = vrotate(force, joint.qoffset) # back to A
    return [szeros(T, 3); force]
    # return spring_relative(:parent, joint, xa, qa, xb, qb; unitary=unitary)
end

@inline function spring_child(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    # spring = unitary ? 1.0 : joint.spring
    # # q = rotation_error(joint, qa, qb, qoff = spring_qoffset(joint))
    # q = spring_extension(joint, xa, qa, xb, qb)
    # distance = spring_distance(joint, q)
    # force = - spring * distance # force in offset frame
    # rotate && (force = vrotate(force, inv(qb) * qa * joint.qoffset)) # rotate back to b frame

    force = spring_force(joint, xa, qa, xb, qb; unitary=unitary)
    # rotate && (force = vrotate(force, inv(qb) * qa * joint.qoffset)) # rotate back to b frame
    rotate && (force = vrotate(force, inv(qb) * qa)) # rotate back to b frame

    return [szeros(T, 3); force]
    # return spring_relative(:child, joint, xa, qa, xb, qb; unitary=unitary)
end
@inline function spring_parent_new(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    return spring_relative(:parent, joint, xa, qa, xb, qb; unitary=unitary)
end

@inline function spring_child_new(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    return spring_relative(:child, joint, xa, qa, xb, qb; unitary=unitary)
end

# vis = Visualizer()
# open(vis)
#
# mech = get_snake(spring=10.0, damper=1.0, Nb=2, gravity=0.0, contact=false)
# initialize!(mech, :snake)
# function ctrl!(mech, k)
#     nu = control_dimension(mech)
#     set_control!(mech, [szeros(6); 11srand(nu-6)]*mech.timestep)
# end
# storage = simulate!(mech, 8.0, ctrl!, record=true, verbose=false)
# visualize(mech, storage, vis=vis)
#
# joint0 = mech.joints[1]
# xa = mech.origin.state.x2[1]
# qa = mech.origin.state.q2[1]
# xb = mech.bodies[1].state.x2[1]
# qb = mech.bodies[1].state.q2[1]
#
#
# spring_parent(joint0.constraints[2], xa, qa, xb, qb)
# spring_child(joint0.constraints[2], xa, qa, xb, qb)
#
# spring_parent_new(joint0.constraints[2], xa, qa, xb, qb)
# spring_child_new(joint0.constraints[2], xa, qa, xb, qb)
#
# impulse_transform(:parent, joint0.constraints[2], xa, qa, xb, qb)
# impulse_transform(:child, joint0.constraints[2], xa, qa, xb, qb)
#

#
# qa = UnitQuaternion(rand(4)...)
# qb = UnitQuaternion(rand(4)...)
# qoff = UnitQuaternion(rand(4)...)
# Vmat(qa \ qb / qoff)
#
# qa \ qb / qoff
# Vmat(orientation_error(joint0.constraints[2], xa, qa, xb, qb))


# #
# # a = 10
# # a = 10
# # a = 10
# # a = 10
# # a = 10

@inline function damper_parent(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector,
        qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where T
    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qoffset)) # rotate back to frame a
    return [szeros(T, 3); force]
end

@inline function damper_child(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector,
        qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where T
    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = - 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qa * qoffset)) # rotate back to frame b
    return [szeros(T, 3); force]
end

spring_parent_jacobian_configuration_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_parent_jacobian_configuration_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuration_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuration_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)

spring_parent_jacobian_velocity_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_parent_jacobian_velocity_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_jacobian_velocity_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_jacobian_velocity_child(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_jacobian_velocity_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_jacobian_velocity_parent(joint::Rotational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)

function spring_parent_jacobian_configuration_parent(joint::Rotational, body1::Node, body2::Node,
        timestep::T; attjac::Bool = true) where T
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    # qoffset = joint.qoffset
    # force = spring_parent(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_parent(joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function spring_parent_jacobian_configuration_child(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    # qoffset = joint.qoffset
    # force = spring_parent(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_parent(joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function spring_child_jacobian_configuration_child(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    # qoffset = joint.qoffset
    # force = spring_child(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_child(joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function spring_child_jacobian_configuration_parent(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωa = current_configuration_velocity(body1.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    # qoffset = joint.qoffset
    # force = spring_child(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_child(joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end


function damper_parent_jacobian_configuration_parent(joint::Rotational, body1::Node, body2::Node,
        timestep::T; attjac::Bool = true) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_parent(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function damper_parent_jacobian_configuration_child(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = current_configuration(body2.state)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_parent(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function damper_child_jacobian_configuration_child(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_child(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat()
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end
function damper_child_jacobian_configuration_parent(joint::Rotational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_child(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb))
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return timestep * [Z; X Q]
end


function spring_parent_jacobian_velocity_parent(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return timestep * [szeros(T, 3, 6); V Ω]
end
function spring_parent_jacobian_velocity_child(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return timestep * [szeros(T, 3, 6); V Ω]
end
function spring_child_jacobian_velocity_child(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return timestep * [szeros(T, 3, 6); V Ω]
end
function spring_child_jacobian_velocity_parent(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return timestep * [szeros(T, 3, 6); V Ω]
end


function damper_parent_jacobian_velocity_parent(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_parent(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return timestep * [szeros(T, 3, 6); V Ω]
end
function damper_parent_jacobian_velocity_child(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_parent(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return timestep * [szeros(T, 3, 6); V Ω]
end
function damper_child_jacobian_velocity_child(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_child(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return timestep * [szeros(T, 3, 6); V Ω]
end
function damper_child_jacobian_velocity_parent(joint::Rotational, body1::Node, body2::Node, timestep::T) where T
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = current_configuration_velocity(body1.state)
    _, _, _, ωb = current_configuration_velocity(body2.state)
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    qoffset = joint.qoffset
    force = damper_child(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return timestep * [szeros(T, 3, 6); V Ω]
end
