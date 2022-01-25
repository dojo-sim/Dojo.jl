## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion; qoff::UnitQuaternion = spring_qoffset(joint))
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
    A = nullspacemat(joint)
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
    # We need to convert the joint.spring_offset which is some minimal coordinate representation of the joint (i.e. axis angle). For us it's the axis angle reresentation.
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    aa = Aᵀ * joint.spring_offset # axis angle
    qoff = axisangle2quaternion(aa)
    return qoff
end

### Spring and damper
## Discrete-time position wrappers (for dynamics)
springforcea(joint::Rotational, bodya::Node, bodyb::Node, Δt; unitary::Bool=false) =
    Δt * springforcea(joint, posargs2(bodya.state)[2], posargs2(bodyb.state)[2], unitary=unitary)
springforceb(joint::Rotational, bodya::Node, bodyb::Node, Δt; unitary::Bool=false) =
    Δt * springforceb(joint, posargs2(bodya.state)[2], posargs2(bodyb.state)[2], unitary=unitary)
damperforcea(joint::Rotational, bodya::Node, bodyb::Node, Δt; unitary::Bool=false) =
    Δt * damperforcea(joint, posargs2(bodya.state)[2], bodya.state.ϕsol[2], posargs2(bodyb.state)[2], bodyb.state.ϕsol[2], unitary=unitary)
damperforceb(joint::Rotational, bodya::Node, bodyb::Node, Δt; unitary::Bool=false) =
    Δt * damperforceb(joint, posargs2(bodya.state)[2], bodya.state.ϕsol[2], posargs2(bodyb.state)[2], bodyb.state.ϕsol[2], unitary=unitary)

springforcea(joint::Rotational3{T}, bodya::Node, bodyb::Node, Δt) where {T} = szeros(T, 6)
springforceb(joint::Rotational3{T}, bodya::Node, bodyb::Node, Δt) where {T} = szeros(T, 6)
damperforcea(joint::Rotational3{T}, bodya::Node, bodyb::Node, Δt) where {T} = szeros(T, 6)
damperforceb(joint::Rotational3{T}, bodya::Node, bodyb::Node, Δt) where {T} = szeros(T, 6)

# Used in energy computation
springforcea(joint::Rotational3{T}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where {T} = szeros(T, 6)
springforceb(joint::Rotational3{T}, qa::UnitQuaternion, qb::UnitQuaternion; rotate::Bool=true, unitary::Bool=false) where {T} = szeros(T, 6)

### Spring and damper
# Force applied by body b on body a expressed in frame a
@inline function springforcea(joint::Rotational{T}, qa::UnitQuaternion, qb::UnitQuaternion;
        rotate::Bool=true, unitary::Bool=false) where {T}
    spring = unitary ? 1.0 : joint.spring
    q = gc(joint, qa, qb, qoff = spring_qoffset(joint))
    distance = spring_distance(joint, q)
    force = spring * distance # force in offset frame
    rotate && (force = vrotate(force, joint.qoffset)) # rotate back to a frame
    return [szeros(T, 3); force]
end

# Force applied by body a on body b expressed in frame b
@inline function springforceb(joint::Rotational{T}, qa::UnitQuaternion, qb::UnitQuaternion;
        rotate::Bool=true, unitary::Bool=false) where {T}
    spring = unitary ? 1.0 : joint.spring
    q = gc(joint, qa, qb, qoff = spring_qoffset(joint))
    distance = spring_distance(joint, q)
    force = - spring * distance # force in offset frame
    rotate && (force = vrotate(force, inv(qb) * qa * joint.qoffset)) # rotate back to b frame
    return [szeros(T, 3); force]
end

# Force applied by body b on body a expressed in frame a
@inline function damperforcea(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector,
        qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where {T}
    damper = unitary ? 1.0 : joint.damper
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qoffset)) # rotate back to frame a
    return [szeros(T, 3); force]
end

# Force applied by body a on body b expressed in frame b
@inline function damperforceb(joint::Rotational{T}, qa::UnitQuaternion, ωa::AbstractVector,
        qb::UnitQuaternion, ωb::AbstractVector; rotate::Bool=true, unitary::Bool=false) where {T}
    damper = unitary ? 1.0 : joint.damper
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    qoffset = joint.qoffset
    velocity = A * (vrotate(ωb, qa \ qb / qoffset) - vrotate(ωa, inv(qoffset))) # in offset frame
    force = - 2 * Aᵀ * A * damper * Aᵀ * velocity # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, inv(qb) * qa * qoffset)) # rotate back to frame b
    return [szeros(T, 3); force]
end


∂springforcea∂posa(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforcea∂posa(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforcea∂posb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforcea∂posb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforceb∂posb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforceb∂posb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂springforceb∂posa(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
∂damperforceb∂posa(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
# ∂springforceb∂posb(joint::Rotational3{T}, body1::Origin, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)
# ∂damperforceb∂posb(joint::Rotational3{T}, body1::Origin, body2::Node, Δt::T; attjac::Bool = true) where T = attjac ? szeros(T, 6, 6) : szeros(T, 6, 7)

∂springforcea∂vela(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂damperforcea∂vela(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂springforcea∂velb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂damperforcea∂velb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂springforceb∂velb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂damperforceb∂velb(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂springforceb∂vela(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
∂damperforceb∂vela(joint::Rotational3{T}, body1::Node, body2::Node, Δt::T) where T = szeros(T, 6, 6)
# ∂springforceb∂velb(joint::Rotational3{T}, body1::Origin, body2::Node, Δt::T) where T = szeros(T, 6, 6)
# ∂damperforceb∂velb(joint::Rotational3{T}, body1::Origin, body2::Node, Δt::T) where T = szeros(T, 6, 6)

function ∂springforcea∂posa(joint::Rotational, body1::Node, body2::Node,
        Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforcea(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> springforcea(joint, UnitQuaternion(qa..., false), qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂damperforcea∂posa(joint::Rotational, body1::Node, body2::Node,
        Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂springforcea∂posb(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforcea(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, qoffset) * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> springforcea(joint, qa, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂damperforcea∂posb(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    x1b, q1b = posargs2(body2.state)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂springforceb∂posb(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforceb(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * Lmat(inv(qa)) * LVᵀmat(qb)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat() * LVᵀmat(qb)
    Q = FiniteDiff.finite_difference_jacobian(qb -> springforceb(joint, qa, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end
function ∂damperforceb∂posb(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(inv(qoffset)) * Lmat(inv(qa))
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qa * qoffset) * Tmat()
    attjac && (Q *= LVᵀmat(qb))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂springforceb∂posa(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # _, _, _, ωa = fullargssol(body1.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    # qoffset = joint.qoffset
    # force = springforceb(joint, qa, qb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    # Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(qb * inv(qoffset)) * Tmat() * LVᵀmat(qa)
    # Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb)) * LVᵀmat(qa)
    Q = FiniteDiff.finite_difference_jacobian(qa -> springforceb(joint, UnitQuaternion(qa..., false), qb)[SVector{3,Int}(4,5,6)], [qa.w, qa.x, qa.y, qa.z])
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

function ∂damperforceb∂posa(joint::Rotational, body1::Node, body2::Node, Δt::T; attjac::Bool = true) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    X = szeros(T, 3, 3)
    Q = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qa \ qb / qoffset) * Rmat(qb * inv(qoffset)) * Tmat()
    Q += ∂vrotate∂q(force, inv(qb) * qa * qoffset) * Rmat(qoffset) * Lmat(inv(qb))
    attjac && (Q *= LVᵀmat(qa))
    Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
    return Δt * [Z; X Q]
end

# function ∂springforceb∂posb(joint::Rotational, body1::Origin, body2::Node, Δt::T; attjac::Bool = true) where T
#     # A = nullspacemat(joint)
#     # Aᵀ = zerodimstaticadjoint(A)
#     # _, _, _, ωb = fullargssol(body2.state)
#     xb, qb = posargs2(body2.state)
#     # qoffset = joint.qoffset
#     # force = springforceb(joint, qb; rotate = false)[SVector{3,Int}(4,5,6)]
#     X = szeros(T, 3, 3)
#     # Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -1.0 * Aᵀ * A * joint.spring * Aᵀ * A * VRmat(inv(qoffset)) * LVᵀmat(qb)
#     # Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat() * LVᵀmat(qb)
#     Q = FiniteDiff.finite_difference_jacobian(qb -> springforceb(joint, UnitQuaternion(qb..., false))[SVector{3,Int}(4,5,6)], [qb.w, qb.x, qb.y, qb.z])
#     attjac && (Q *= LVᵀmat(qb))
#     Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
#     return Δt * [Z; X Q]
# end
# function ∂damperforceb∂posb(joint::Rotational, body1::Origin, body2::Node, Δt::T; attjac::Bool = true) where T
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     _, _, _, ωb = fullargssol(body2.state)
#     xb, qb = posargs2(body2.state)
#     qoffset = joint.qoffset
#     force = damperforceb(joint, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
#     X = szeros(T, 3, 3)
#     Q = ∂vrotate∂p(force, inv(qb) * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂q(ωb, qb / qoffset) * Rmat(inv(qoffset))
#     Q += ∂vrotate∂q(force, inv(qb) * qoffset) * Rmat(qoffset) * Tmat()
#     attjac && (Q *= LVᵀmat(qb))
#     Z = attjac ? szeros(T, 3, 6) : szeros(T, 3, 7)
#     return Δt * [Z; X Q]
# end


function ∂springforcea∂vela(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforcea∂vela(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforcea∂velb(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforcea∂velb(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforcea(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, qoffset) * 2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforceb∂velb(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforceb∂velb(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qa \ qb / qoffset)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂springforceb∂vela(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    V = szeros(T, 3, 3)
    Ω = szeros(T, 3, 3)
    return Δt * [szeros(T, 3, 6); V Ω]
end
function ∂damperforceb∂vela(joint::Rotational, body1::Node, body2::Node, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    _, _, _, ωa = fullargssol(body1.state)
    _, _, _, ωb = fullargssol(body2.state)
    xa, qa = posargs2(body1.state)
    xb, qb = posargs2(body2.state)
    qoffset = joint.qoffset
    force = damperforceb(joint, qa, ωa, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
    V = szeros(T, 3, 3)
    Ω = ∂vrotate∂p(force, inv(qb) * qa * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * -1.0 * ∂vrotate∂p(ωa, inv(qoffset))
    return Δt * [szeros(T, 3, 6); V Ω]
end
# function ∂springforceb∂velb(joint::Rotational, body1::Origin, body2::Node, Δt::T) where T
#     V = szeros(T, 3, 3)
#     Ω = szeros(T, 3, 3)
#     return Δt * [szeros(T, 3, 6); V Ω]
# end
# function ∂damperforceb∂velb(joint::Rotational, body1::Origin, body2::Node, Δt::T) where T
#     A = nullspacemat(joint)
#     Aᵀ = zerodimstaticadjoint(A)
#     _, _, _, ωb = fullargssol(body2.state)
#     xb, qb = posargs2(body2.state)
#     qoffset = joint.qoffset
#     force = damperforceb(joint, qb, ωb; rotate = false)[SVector{3,Int}(4,5,6)]
#     V = szeros(T, 3, 3)
#     Ω = ∂vrotate∂p(force, inv(qb) * qoffset) * -2 * Aᵀ * A * joint.damper * Aᵀ * A * ∂vrotate∂p(ωb, qb / qoffset)
#     return Δt * [szeros(T, 3, 6); V Ω]
# end
