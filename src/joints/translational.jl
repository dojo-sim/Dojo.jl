mutable struct Translational{T,N} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    spring::T
    damper::T
    Fτ::SVector{3,T}
    function Translational{T,N}(body1::AbstractBody, body2::AbstractBody;
            p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3), spring = zero(T), damper = zero(T)
        ) where {T,N}
        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]
        Fτ = zeros(T,3)
        new{T,N}(V3, V12, vertices, spring, damper, Fτ), body1.id, body2.id
    end
end
Translational0 = Translational{T,0} where T
Translational1 = Translational{T,1} where T
Translational2 = Translational{T,2} where T
Translational3 = Translational{T,3} where T
springforcea(joint::Translational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
springforceb(joint::Translational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
springforceb(joint::Translational{T,3}, body1::Origin, body2::Body, Δt::T, childid) where T = szeros(T, 6)
damperforcea(joint::Translational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
damperforceb(joint::Translational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
damperforceb(joint::Translational{T,3}, body1::Origin, body2::Body, Δt::T, childid) where T = szeros(T, 6)
function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Translational{T,N}) where {T,N}
    summary(io, constraint)
    println(io,"")
    println(io, " V3:       "*string(constraint.V3))
    println(io, " V12:      "*string(constraint.V12))
    println(io, " vertices: "*string(constraint.vertices))
end
### Constraints and derivatives
## Position level constraints (for dynamics)
@inline function g(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end
@inline function g(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    point2 = xb + vrotate(joint.vertices[2], qb)

    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * (Lmat(UnitQuaternion(point2)) - Lmat(UnitQuaternion(xa)))

    return X, Q
end
@inline function ∂g∂posb(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end
@inline function ∂g∂posb(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
    X = I
    Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end


# # Derivatives accounting for quaternion specialness
# @inline function ∂g∂ʳposa(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
#     X, Q = ∂g∂posa(joint, xa, qa, xb, qb)
#     Q = Q * LVᵀmat(qa)
#
#     return [X Q]
# end
# @inline function ∂g∂ʳposb(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
#     X, Q = ∂g∂posb(joint, xa, qa, xb, qb)
#     Q = Q * LVᵀmat(qb)
#
#     return [X Q]
# end
# @inline function ∂g∂ʳposb(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
#     X, Q = ∂g∂posb(joint, xb, qb)
#     Q = Q * LVᵀmat(qb)
#
#     return [X Q]
# end

# @inline function g(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
#     vertices = joint.vertices
#     return xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
# end
# @inline function g(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
#     vertices = joint.vertices
#     return xb + vrotate(vertices[2], qb) - vertices[1]
# end
# ## Derivatives NOT accounting for quaternion specialness
# @inline function ∂g∂posa(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
#     X = -I(3)
#     Q = -∂vrotate∂q(joint.vertices[1], qa)
#     return X, Q
# end
# @inline function ∂g∂posb(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
#     X = I(3)
#     Q = ∂vrotate∂q(joint.vertices[2], qb)
#     return X, Q
# end
# @inline function ∂g∂posb(joint::Translational, xb::AbstractVector, qb::UnitQuaternion)
#     X = I(3)
#     Q = ∂vrotate∂q(joint.vertices[2], qb)
#     return X, Q
# end
function ∂g∂ʳposa(joint::Translational{T}, statea::State, stateb::State, Δt) where T
    vertices = joint.vertices
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    gw = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    # ga = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))

    # X = -I(3)
    # Q = rotation_matrix(qa) * (∂vrotate∂p(gw, inv(qa)) * -1.0 * ∂vrotate∂q(vertices[1], qa) + ∂vrotate∂q(gw, inv(qa)) * Tmat()) * LVᵀmat(qa)
    # Q = -rotation_matrix(inv(qa)) * ∂vrotate∂q(joint.vertices[1], qa) * LVᵀmat(qa)
    # Q = transpose(-1.0 * rotation_matrix(inv(qa)) * skew(rotation_matrix(qa) * joint.vertices[1]))
    X = - rotation_matrix(inv(qa))
    Q = skew(vertices[1])
    # @show scn.(Q0 - Q)
    return [X Q]
end
function ∂g∂ʳposb(joint::Translational{T}, statea::State, stateb::State, Δt) where T
    # vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qb))
    vertices = joint.vertices
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    gw = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    # ga = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))

    # X = I(3)
    # Q = rotation_matrix(qa) * ∂vrotate∂p(gw, inv(qa)) * ∂vrotate∂q(vertices[2], qb) * LVᵀmat(qb)
    # Q = rotation_matrix(inv(qa)) * ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    # Q = transpose(rotation_matrix(inv(qb)) * skew(rotation_matrix(qb) * joint.vertices[2]))
    # @show scn.(Q0 - Q)
    X = - rotation_matrix(inv(qa)) * rotation_matrix(qb) * skew(vertices[2]apu)
    Q = skew(vertices[1])
    return [X Q]
end
function ∂g∂ʳposb(joint::Translational{T}, stateb::State, Δt) where T
    # vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qb))
    vertices = joint.vertices
    xb, qb = posargsk(stateb)
    # gw = xb + vrotate(vertices[2], qb) - vertices[1]

    X = I(3)
    Q = ∂vrotate∂q(vertices[2], qb) * LVᵀmat(qb)
    # Q = ∂vrotate∂q(joint.vertices[2], qb) * LVᵀmat(qb)
    # Q = transpose(rotation_matrix(inv(qb)) * skew(rotation_matrix(qb) * joint.vertices[2]))
    # @show scn.(Q0 - Q)
    return [X Q]
end
## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
@inline function ∂2g∂posaa(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
    Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    XX = szeros(T, 9, 3)
    XQ = -kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) - kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = -kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Lpos)*∂L∂qsplit(T) + kron(VLᵀmat(qa)*Ltpos,2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(VLᵀmat(qa),2*VLᵀmat(qa)*Lmat(qb)*Lmat(UnitQuaternion(joint.vertices[2])))*∂Lᵀ∂qsplit(T) + kron(VLᵀmat(qa)*Lmat(qb)*Lᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) + kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = szeros(T, 9, 3)
    QQ = kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb),2*VLᵀmat(qa))*∂R∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb)*Rᵀmat(qa),2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Rmat(qa)*Rᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa)*Rmat(qa))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Translational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VRᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*Vmat(T))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end


### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Translational{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
    F = joint.Fτ
    vertices = joint.vertices
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)

    Fa = vrotate(-F, qa)
    Fb = -Fa

    τa = vrotate(torqueFromForce(Fa, vrotate(vertices[1], qa)),inv(qa)) # in local coordinates
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    statea.Fk[end] += Fa
    statea.τk[end] += τa
    stateb.Fk[end] += Fb
    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end
@inline function applyFτ!(joint::Translational{T}, stateb::State, Δt::T, clear::Bool) where T
    F = joint.Fτ
    vertices = joint.vertices
    _, qb = posargsk(stateb)

    Fb = F
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    stateb.Fk[end] += Fb
    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

## Forcing derivatives (for linearization)
# Control derivatives
@inline function ∂Fτ∂ua(joint::Translational, statea::State, stateb::State, Δt::T) where T
    vertices = joint.vertices
    _, qa = posargsk(statea)

    BFa = -VLmat(qa) * RᵀVᵀmat(qa)
    Bτa = -skew(vertices[1])

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Translational, statea::State, stateb::State, Δt::T) where T
    vertices = joint.vertices
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    qbinvqa = qb\qa

    BFb = VLmat(qa) * RᵀVᵀmat(qa)
    Bτb = skew(vertices[2]) * VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Translational, stateb::State, Δt::T) where T
    vertices = joint.vertices
    _, qb = posargsk(stateb)

    BFb = I
    Bτb = skew(vertices[2]) * VLᵀmat(qb) * RVᵀmat(qb)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Translational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXa = szeros(T,3,3)
    FaQa = -2*VRᵀmat(qa)*Rmat(UnitQuaternion(F))
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,4)
    FbXa = szeros(T,3,3)
    FbQa = 2*VRᵀmat(qa)*Rmat(UnitQuaternion(F))
    τbXa = szeros(T,3,3)
    τbQa = 2*skew(vertices[2])*VLᵀmat(qb)*Rmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(F))

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Translational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,4)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,4)
    τbXb = szeros(T,3,3)
    τbQb = 2*skew(vertices[2])*VLᵀmat(qb)*Lmat(qa)*Lmat(UnitQuaternion(F))*Lᵀmat(qa)#*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Translational{T}, stateb::State, Δt::T) where T
    xb, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,4)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,4)
    τbXb = szeros(T,3,3)
    τbQb = 2*skew(vertices[2])*VLᵀmat(qb)*Lmat(UnitQuaternion(F))#*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end


### Minimal coordinates
## Position and velocity offsets
@inline function getPositionDelta(joint::Translational, body1::AbstractBody, body2::Body, x::SVector)
    Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::Translational, body1::AbstractBody, body2::Body, v::SVector)
    Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
    return Δv
end

## Minimal coordinate calculation
@inline function minimalCoordinates(joint::Translational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::Translational, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
end
@inline function minimalVelocities(joint::Translational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * (stateb.vc - statea.vc)
end
@inline function minimalVelocities(joint::Translational, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * stateb.vc
end
