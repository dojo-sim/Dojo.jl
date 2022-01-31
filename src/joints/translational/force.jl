spring_parent(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_parent(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
spring_child(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * spring_child(joint, current_configuration(bodya.state)..., current_configuration(bodyb.state)..., unitary=unitary)
damper_parent(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_parent(joint, current_configuration(bodya.state)..., bodya.state.vsol[2],
    bodya.state.ϕsol[2], current_configuration(bodyb.state)..., bodyb.state.vsol[2], bodyb.state.ϕsol[2], unitary=unitary)
damper_child(joint::Translational, bodya::Node, bodyb::Node, timestep; unitary::Bool=false) =
    timestep * damper_child(joint, current_configuration(bodya.state)..., bodya.state.vsol[2],
    bodya.state.ϕsol[2], current_configuration(bodyb.state)..., bodyb.state.vsol[2], bodyb.state.ϕsol[2], unitary=unitary)

spring_parent(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
spring_child(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_parent(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)
damper_child(joint::Translational3{T}, bodya::Node, bodyb::Node, timestep) where T = szeros(T, 6)

spring_parent(joint::Translational3{T}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)
spring_child(joint::Translational3{T}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T = szeros(T, 6)

@inline function spring_parent(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    # spring = unitary ? 1.0 : joint.spring
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # distance = A * position_error(joint, xa, qa, xb, qb) .- joint.spring_offset
    # force = spring * Aᵀ * distance # Currently assumes same spring constant in all directions
    # forceA = force # in the A frame
    # rotate && (force = vrotate(force, qa)) # rotate back to world frame
    # torque = skew(joint.vertices[1]) * forceA

    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δmincoord0 = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb) # in the a frame
    forceA0 = spring * Aᵀ * Δmincoord0 # in the a frame
    Fτ0 = impulse_transform_parent(joint, xa, qa, xb, qb) * forceA0
    # @show norm(Fτ0 - [force; torque], Inf)

    # return [force; torque]
    return Fτ0
end

@inline function spring_child(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where T
    spring = unitary ? 1.0 : joint.spring
    # A = nullspace_mask(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # distance = A * position_error(joint, xa, qa, xb, qb) .- joint.spring_offset
    # force = - spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    # forceA = force
    # rotate && (force = vrotate(force, qa)) # rotate back to world frame
    #
    # pa_a = rotation_matrix(inv(qa)) * (xa + rotation_matrix(qa) * joint.vertices[1]) # body a kinematics point
    # cb_a = rotation_matrix(inv(qa)) * (xb) # body b com
    # ra = pa_a - cb_a
    # torque = rotation_matrix(inv(qb) * qa) * skew(ra) * forceA

    spring = unitary ? 1.0 : joint.spring
    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    Δmincoord0 = joint.spring_offset .- minimal_coordinates(joint, xa, qa, xb, qb) # in the a frame
    forceA0 = spring * Aᵀ * Δmincoord0 # in the a frame
    Fτ0 = impulse_transform_child(joint, xa, qa, xb, qb) * forceA0
    # @show norm(Fτ0 - [force; torque], Inf)

    # return [force; torque]
    return Fτ0
end

@inline function damper_parent(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where T
    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)

    pa_b = rotation_matrix(inv(qb)) * (xa + rotation_matrix(qa) * joint.vertices[1]) # body a kinematics point
    cb_b = xb # body b com
    rb = pa_b - cb_b
    vpb = vb + vrotate(skew(ωb) * rb, qb)
    vpa = va + vrotate(skew(ωa) * joint.vertices[1], qa)

    # velocity = A * vrotate(vpb - vpa, inv(qa))
    velocity = A * vrotate(vb - va, inv(qa))
    force = damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    forceA = force # in the A frame
    rotate && (force = vrotate(force, qa)) # rotate back to world frame

    torque = skew(joint.vertices[1]) * forceA




    return [force; torque]
end

@inline function damper_child(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        xb::AbstractVector, qb::UnitQuaternion, vb::AbstractVector, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where T
    damper = unitary ? 1.0 : joint.damper
    A = nullspace_mask(joint)
    Aᵀ = zerodimstaticadjoint(A)

    pa_b = rotation_matrix(inv(qb)) * (xa + rotation_matrix(qa) * joint.vertices[1]) # body a kinematics point
    cb_b = xb # body b com
    rb = pa_b - cb_b
    vpb = vb + vrotate(skew(ωb) * rb, qb)
    vpa = va + vrotate(skew(ωa) * joint.vertices[1], qa)

    # velocity = A * vrotate(vpb - vpa, inv(qa))
    velocity = A * vrotate(vb - va, inv(qa))
    force = - damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    forceA = force
    rotate && (force = vrotate(force, qa)) # rotate back to world frame

    pa_a = rotation_matrix(inv(qa)) * (xa + rotation_matrix(qa) * joint.vertices[1]) # body a kinematics point
    cb_a = rotation_matrix(inv(qa)) * (xb) # body b com
    ra = pa_a - cb_a
    torque = rotation_matrix(inv(qb) * qa) * skew(ra) * forceA
    return [force; torque]
end

spring_parent_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_parent_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
spring_child_jacobian_configuraion_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_parent_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
damper_child_jacobian_configuration_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)

spring_parent_jacobian_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_parent_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
spring_child_configuration_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
# spring_child_jacobian_velocity_child(joint::Translational3{T}, body1::Origin, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_parent_jacobian_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_configuration_velocity_child(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
damper_child_configuration_velocity_parent(joint::Translational3{T}, body1::Node, body2::Node, timestep::T) where T = szeros(T, 6, 6)
# damper_child_configuration_velocity_child(joint::Translational3{T}, body1::Origin, body2::Node, timestep::T) where T = szeros(T, 6, 6)

function spring_parent_jacobian_configuration_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)

    X = FiniteDiff.finite_difference_jacobian(xa -> spring_parent(joint, xa, qa, xb, qb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_parent(joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end

function damper_parent_jacobian_configuration_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]

    X = FiniteDiff.finite_difference_jacobian(xa -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> damper_parent(joint, xa, UnitQuaternion(qa..., false), va, ωa, xb, qb, vb, ωb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end

function spring_parent_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)

    X = FiniteDiff.finite_difference_jacobian(xb -> spring_parent(joint, xa, qa, xb, qb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_parent(joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end

function damper_parent_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]

    X = FiniteDiff.finite_difference_jacobian(xb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> damper_parent(joint, xa, qa, va, ωa, xb, UnitQuaternion(qb..., false), vb, ωb), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end

function spring_child_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)

    X = FiniteDiff.finite_difference_jacobian(xb -> spring_child(joint, xa, qa, xb, qb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> spring_child(joint, xa, qa, xb, UnitQuaternion(qb..., false)), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end

function damper_child_jacobian_configuration_child(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]

    X = FiniteDiff.finite_difference_jacobian(xb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> damper_child(joint, xa, qa, va, ωa, xb, UnitQuaternion(qb..., false), vb, ωb), [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    return timestep * [X Q]
end

function spring_child_jacobian_configuraion_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)

    X = FiniteDiff.finite_difference_jacobian(xa -> spring_child(joint, xa, qa, xb, qb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> spring_child(joint, xa, UnitQuaternion(qa..., false), xb, qb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end

function damper_child_jacobian_configuration_parent(joint::Translational, body1::Node, body2::Node, timestep::T; attjac::Bool = true) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    ωb = body2.state.ϕsol[2]

    X = FiniteDiff.finite_difference_jacobian(xa -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), xa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> damper_child(joint, xa, UnitQuaternion(qa..., false), va, ωa, xb, qb, vb, ωb), [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    return timestep * [X Q]
end






function spring_parent_jacobian_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return timestep * szeros(T, 6, 6)
end

function spring_parent_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return timestep * szeros(T, 6, 6)
end

function spring_child_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return timestep * szeros(T, 6, 6)
end

function spring_child_configuration_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    return timestep * szeros(T, 6, 6)
end


function damper_parent_jacobian_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]

    V = FiniteDiff.finite_difference_jacobian(va -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), va)
    Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωa)
    return timestep * [V Ω]
end

function damper_parent_jacobian_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]

    V = FiniteDiff.finite_difference_jacobian(vb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), vb)
    Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_parent(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωb)
    return timestep * [V Ω]
end

function damper_child_configuration_velocity_child(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]

    V = FiniteDiff.finite_difference_jacobian(vb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), vb)
    Ω = FiniteDiff.finite_difference_jacobian(ωb -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωb)
    return timestep * [V Ω]
end

function damper_child_configuration_velocity_parent(joint::Translational, body1::Node, body2::Node, timestep::T) where T
    xa, qa = current_configuration(body1.state)
    xb, qb = current_configuration(body2.state)
    va = body1.state.vsol[2]
    ωa = body1.state.ϕsol[2]
    vb = body2.state.vsol[2]
    ωb = body2.state.ϕsol[2]

    V = FiniteDiff.finite_difference_jacobian(va -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), va)
    Ω = FiniteDiff.finite_difference_jacobian(ωa -> damper_child(joint, xa, qa, va, ωa, xb, qb, vb, ωb), ωa)
    return timestep * [V Ω]
end
