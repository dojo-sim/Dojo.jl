## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    rotate && (d = vrotate(d, inv(qa)))
    return d
end
@inline function gc(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end

### Spring and damper
## Forces for dynamics
## Discrete-time position wrappers (for dynamics)
springforcea(joint::Translational, statea::State, stateb::State, Δt) = Δt * springforcea(joint, posargsk(statea)..., posargsk(stateb)...)
springforceb(joint::Translational, statea::State, stateb::State, Δt) = Δt * springforceb(joint, posargsk(statea)..., posargsk(stateb)...)
springforceb(joint::Translational, stateb::State, Δt) = Δt * springforceb(joint, posargsk(stateb)...)
damperforcea(joint::Translational, statea::State, stateb::State, Δt) = Δt * damperforcea(joint, posargsk(statea)[2], statea.vsol[2], stateb.vsol[2])
damperforceb(joint::Translational, statea::State, stateb::State, Δt) = Δt * damperforceb(joint, posargsk(statea)[2], statea.vsol[2], stateb.vsol[2])
damperforceb(joint::Translational, stateb::State, Δt) = Δt * damperforceb(joint, stateb.vsol[2])

springforcea(joint::Translational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
springforceb(joint::Translational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
springforceb(joint::Translational{T,3}, stateb::State, Δt) where {T} = szeros(T, 6)
damperforcea(joint::Translational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
damperforceb(joint::Translational{T,3}, statea::State, stateb::State, Δt) where {T} = szeros(T, 6)
damperforceb(joint::Translational{T,3}, stateb::State, Δt) where {T} = szeros(T, 6)


### Spring and damper
## Forces for dynamics
# Force applied by body b on body a expressed in world frame
@inline function springforcea(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = joint.spring * Aᵀ * distance # Currently assumes same spring constant in all directions
    rotate && (force = vrotate(force, qa)) # rotate back to world frame
    # return 0.0 * [force; szeros(T, 3)]
    return [force; szeros(T, 3)]
end
# Force applied by body a on body b expressed in world frame
@inline function springforceb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = - joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    rotate && (force = vrotate(force, qa)) # rotate back to world frame
    # return 0.0 * [force; szeros(T, 3)]
    return [force; szeros(T, 3)]
end
# Force applied by origin on body b expressed in world frame
@inline function springforceb(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return 0.0 * [force; szeros(T, 3)]
end

# Force applied by body b on body a expressed in world frame
@inline function damperforcea(joint::Translational{T}, qa::UnitQuaternion, va::AbstractVector, vb::AbstractVector; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vrotate(vb - va, inv(qa))
    force = Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qa)) # rotate back to world frame
    return 0.0 * [force; szeros(T, 3)]
end
# Force applied by body a on body b expressed in world frame
@inline function damperforceb(joint::Translational{T}, qa::UnitQuaternion, va::AbstractVector, vb::AbstractVector; rotate::Bool = true) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vrotate(vb - va, inv(qa))
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    rotate && (force = vrotate(force, qa)) # rotate back to world frame
    return 0.0 * [force; szeros(T, 3)]
end
# Force applied by origin on body b expressed in world frame
@inline function damperforceb(joint::Translational{T}, vb::AbstractVector) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vb
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return 0.0 * [force; szeros(T, 3)]
end

∂springforcea∂posa(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforcea∂posa(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforcea∂posb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforcea∂posb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂posb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂posb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂posa(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂posa(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂posb(joint::Translational{T,3}, body1::Origin, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂posb(joint::Translational{T,3}, body1::Origin, body2::Body, childid) where T = szeros(T, 6, 6)

∂springforcea∂vela(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforcea∂vela(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforcea∂velb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforcea∂velb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂velb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂velb(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂vela(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂vela(joint::Translational{T,3}, body1::Body, body2::Body, childid) where T = szeros(T, 6, 6)
∂springforceb∂velb(joint::Translational{T,3}, body1::Origin, body2::Body, childid) where T = szeros(T, 6, 6)
∂damperforceb∂velb(joint::Translational{T,3}, body1::Origin, body2::Body, childid) where T = szeros(T, 6, 6)

function ∂springforcea∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)

    X = ForwardDiff.jacobian(xa -> springforcea(joint, xa, qa, xb, qb)[SVector{3,Int}(1,2,3)], xa)
    Q = ForwardDiff.jacobian(qa -> springforcea(joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(1,2,3)], [qa.w, qa.x, qa.y, qa.z]) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforcea∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    X = szeros(T, 3, 3)
    Q = ForwardDiff.jacobian(qa -> damperforcea(joint, UnitQuaternion(qa..., false), va, vb)[SVector{3,Int}(1,2,3)], [qa.w, qa.x, qa.y, qa.z]) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforcea∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # xb, qb = posargsk(body2.state)
    # X = Aᵀ * A * joint.spring * Aᵀ * A
    # Q = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)

    X = ForwardDiff.jacobian(xb -> springforcea(joint, xa, qa, xb, qb)[SVector{3,Int}(1,2,3)], xb)
    Q = ForwardDiff.jacobian(qb -> springforcea(joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(1,2,3)], [qb.w, qb.x, qb.y, qb.z]) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforcea∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # xb, qb = posargsk(body2.state)
    # X = - Aᵀ * A * joint.spring * Aᵀ * A
    # Q = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)

    X = ForwardDiff.jacobian(xb -> springforceb(joint, xa, qa, xb, qb)[SVector{3,Int}(1,2,3)], xb)
    Q = ForwardDiff.jacobian(qb -> springforceb(joint, xa, qa, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(1,2,3)], [qb.w, qb.x, qb.y, qb.z]) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # xa, qa = posargsk(body1.state)
    # X = Aᵀ * A * joint.spring * Aᵀ * A
    # Q = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[1], qa) * LVᵀmat(qa)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)

    X = ForwardDiff.jacobian(xa -> springforceb(joint, xa, qa, xb, qb)[SVector{3,Int}(1,2,3)], xa)
    Q = ForwardDiff.jacobian(qa -> springforceb(joint, xa, UnitQuaternion(qa..., false), xb, qb)[SVector{3,Int}(1,2,3)], [qa.w, qa.x, qa.y, qa.z]) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # X = szeros(T, 3, 3)
    # Q = szeros(T, 3, 3)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    X = szeros(T, 3, 3)
    Q = ForwardDiff.jacobian(qa -> damperforceb(joint, UnitQuaternion(qa..., false), va, vb)[SVector{3,Int}(1,2,3)], [qa.w, qa.x, qa.y, qa.z]) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # xb, qb = posargsk(body2.state)
    # X = - Aᵀ * A * joint.spring * Aᵀ * A
    # Q = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    xb, qb = posargsk(body2.state)

    X = ForwardDiff.jacobian(xb -> springforceb(joint, xb, qb)[SVector{3,Int}(1,2,3)], xb)
    Q = ForwardDiff.jacobian(qb -> springforceb(joint, xb, UnitQuaternion(qb..., false))[SVector{3,Int}(1,2,3)], [qb.w, qb.x, qb.y, qb.z]) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end

function ∂springforcea∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    return Δt * szeros(T, 6, 6)
end
function ∂damperforcea∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # V = - Aᵀ * A * joint.damper * Aᵀ * A
    # Ω = szeros(T, 3, 3)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    V = ForwardDiff.jacobian(va -> damperforcea(joint, qa, va, vb)[SVector{3,Int}(1,2,3)], va)
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforcea∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    return Δt * szeros(T, 6, 6)
end
function ∂damperforcea∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # V = Aᵀ * A * joint.damper * Aᵀ * A
    # Ω = szeros(T, 3, 3)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    V = ForwardDiff.jacobian(vb -> damperforcea(joint, qa, va, vb)[SVector{3,Int}(1,2,3)], vb)
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # V = - Aᵀ * A * joint.damper * Aᵀ * A
    # Ω = szeros(T, 3, 3)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    V = ForwardDiff.jacobian(vb -> damperforceb(joint, qa, va, vb)[SVector{3,Int}(1,2,3)], vb)
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # V = Aᵀ * A * joint.damper * Aᵀ * A
    # Ω = szeros(T, 3, 3)
    xa, qa = posargsk(body1.state)
    xb, qb = posargsk(body2.state)
    va = body1.state.vsol[2]
    vb = body2.state.vsol[2]

    V = ForwardDiff.jacobian(va -> damperforceb(joint, qa, va, vb)[SVector{3,Int}(1,2,3)], va)
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂velb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂velb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # V = - Aᵀ * A * joint.damper * Aᵀ * A
    # Ω = szeros(T, 3, 3)
    xb, qb = posargsk(body2.state)
    vb = body2.state.vsol[2]

    V = ForwardDiff.jacobian(vb -> damperforceb(joint, vb)[SVector{3,Int}(1,2,3)], vb)
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
