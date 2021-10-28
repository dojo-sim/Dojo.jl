## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
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
damperforcea(joint::Translational, statea::State, stateb::State, Δt) = Δt * damperforcea(joint, statea.vsol[2], stateb.vsol[2])
damperforceb(joint::Translational, statea::State, stateb::State, Δt) = Δt * damperforceb(joint, statea.vsol[2], stateb.vsol[2])
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
    xb::AbstractVector, qb::UnitQuaternion) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return [force; szeros(T, 3)]
end
# Force applied by body a on body b expressed in world frame
@inline function springforceb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return [force; szeros(T, 3)]
end
# Force applied by origin on body b expressed in world frame
@inline function springforceb(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return [force; szeros(T, 3)]
end

# Force applied by body b on body a expressed in world frame
@inline function damperforcea(joint::Translational{T}, va::AbstractVector, vb::AbstractVector) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * (vb - va)
    force = Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return [force; szeros(T, 3)]
end
# Force applied by body a on body b expressed in world frame
@inline function damperforceb(joint::Translational{T}, va::AbstractVector, vb::AbstractVector) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * (vb - va)
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return [force; szeros(T, 3)]
end
# Force applied by origin on body b expressed in world frame
@inline function damperforceb(joint::Translational{T}, vb::AbstractVector) where {T}
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vb
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return [force; szeros(T, 3)]
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
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, qa = posargsk(body1.state)
    X = - Aᵀ * A * joint.spring * Aᵀ * A
    Q = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[1], qa) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforcea∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforcea∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xb, qb = posargsk(body2.state)
    X = Aᵀ * A * joint.spring * Aᵀ * A
    Q = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforcea∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xb, qb = posargsk(body2.state)
    X = - Aᵀ * A * joint.spring * Aᵀ * A
    Q = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xa, qa = posargsk(body1.state)
    X = Aᵀ * A * joint.spring * Aᵀ * A
    Q = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[1], qa) * LVᵀmat(qa)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posa(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂springforceb∂posb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    xb, qb = posargsk(body2.state)
    X = - Aᵀ * A * joint.spring * Aᵀ * A
    Q = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    return Δt * [X Q; szeros(T, 3, 6)]
end
function ∂damperforceb∂posb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    X = szeros(T, 3, 3)
    Q = szeros(T, 3, 3)
    return Δt * [X Q; szeros(T, 3, 6)]
end

function ∂springforcea∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # x1a, q1a = posargsk(body1.state)
    # _, _, _, ωa = fullargssol(body1.state)
    # # xa, qa = posargsnext(body1.state, Δt)
    # xa, qa = posargshalf(body1.state, Δt)
    # # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt
    # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt/2
    # Ω = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[1], qa) * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    # # return Δt * [V Ω; szeros(T, 3, 6)]
    return Δt * szeros(T, 6, 6)
end
function ∂damperforcea∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforcea∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # x1b, q1b = posargsk(body2.state)
    # _, _, _, ωb = fullargssol(body2.state)
    # # xb, qb = posargsnext(body2.state, Δt)
    # xb, qb = posargshalf(body2.state, Δt)
    # # V = Aᵀ * A * joint.spring * Aᵀ * A * Δt
    # V = Aᵀ * A * joint.spring * Aᵀ * A * Δt/2
    # Ω = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    # # return Δt * [V Ω; szeros(T, 3, 6)]
    return Δt * szeros(T, 6, 6)
end
function ∂damperforcea∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # x1b, q1b = posargsk(body2.state)
    # _, _, _, ωb = fullargssol(body2.state)
    # # xb, qb = posargsnext(body2.state, Δt)
    # xb, qb = posargshalf(body2.state, Δt)
    # # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt
    # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt/2
    # Ω = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    # # return Δt * [V Ω; szeros(T, 3, 6)]
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂velb(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # x1a, q1a = posargsk(body1.state)
    # _, _, _, ωa = fullargssol(body1.state)
    # # xa, qa = posargsnext(body1.state, Δt)
    # xa, qa = posargshalf(body1.state, Δt)
    # # V = Aᵀ * A * joint.spring * Aᵀ * A * Δt
    # V = Aᵀ * A * joint.spring * Aᵀ * A * Δt/2
    # Ω = Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[1], qa) * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    # # return Δt * [V Ω; szeros(T, 3, 6)]
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂vela(joint::Translational, body1::Body, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
function ∂springforceb∂velb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    # A = nullspacemat(joint)
    # Aᵀ = zerodimstaticadjoint(A)
    # x1b, q1b = posargsk(body2.state)
    # _, _, _, ωb = fullargssol(body2.state)
    # # xb, qb = posargsnext(body2.state, Δt)
    # xb, qb = posargshalf(body2.state, Δt)
    # # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt
    # V = - Aᵀ * A * joint.spring * Aᵀ * A * Δt/2
    # Ω = - Aᵀ * A * joint.spring * Aᵀ * A * ∂vrotate∂q(joint.vertices[2], qb) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    # # return Δt * [V Ω; szeros(T, 3, 6)]
    return Δt * szeros(T, 6, 6)
end
function ∂damperforceb∂velb(joint::Translational, body1::Origin, body2::Body, Δt::T) where T
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    V = - Aᵀ * A * joint.damper * Aᵀ * A
    Ω = szeros(T, 3, 3)
    return Δt * [V Ω; szeros(T, 3, 6)]
end
