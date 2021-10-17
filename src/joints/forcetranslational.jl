mutable struct ForceTranslational12{T,N} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames

    spring::T
    damper::T

    Fτ::SVector{3,T}

    function ForceTranslational12{T,N}(body1::AbstractBody, body2::AbstractBody;
            p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3), spring = zero(T), damper = zero(T)
        ) where {T,N}

        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)

        new{T,N}(V3, V12, vertices, spring, damper, Fτ), body1.id, body2.id
    end
end

ForceTranslational120 = ForceTranslational12{T,0} where T
ForceTranslational121 = ForceTranslational12{T,1} where T
ForceTranslational122 = ForceTranslational12{T,2} where T
ForceTranslational123 = ForceTranslational12{T,3} where T

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::ForceTranslational12{T,N}) where {T,N}
    summary(io, constraint)
    println(io,"")
    println(io, " V3:       "*string(constraint.V3))
    println(io, " V12:      "*string(constraint.V12))
    println(io, " vertices: "*string(constraint.vertices))
end


### Constraints and derivatives
## Position level constraint wrappers
@inline function g(joint::ForceTranslational12, body1::Body, body2::Body, Δt) #USELESS
    Ac = constraintmat(joint)
    An = nullspacemat(joint)
    Acᵀ = zerodimstaticadjoint(Ac)
    Anᵀ = zerodimstaticadjoint(An)
    Fτ = springforce(joint, body1.state, body2.state, Δt) + damperforce(joint, body1.state, body2.state, Δt)
    return Acᵀ * Ac * gc(joint, body1.state, body2.state, Δt) + Anᵀ * An * Fτ
end

@inline function g(joint::ForceTranslational12, body1::Origin, body2::Body, Δt) #USELESS
    Ac = constraintmat(joint)
    An = nullspacemat(joint)
    Acᵀ = zerodimstaticadjoint(Ac)
    Anᵀ = zerodimstaticadjoint(An)
    Fτ = springforce(joint, body2.state, Δt) + damperforce(joint, body2.state, Δt)
    return Acᵀ * Ac * gc(joint, body2.state, Δt) + Anᵀ * An * Fτ
end

@inline function g(joint::ForceTranslational12, body1::Body, body2::Body)
    Ac = constraintmat(joint)
    An = nullspacemat(joint)
    Acᵀ = zerodimstaticadjoint(Ac)
    Anᵀ = zerodimstaticadjoint(An)
    Fτ = springforce(joint, body1.state, body1.state) + damperforce(joint, body1.state, body2.state)
    @show "here"
    return Acᵀ * Ac * gc(joint, body1.state, body2.state) + Anᵀ * An * Fτ
end

@inline function g(joint::ForceTranslational12, body1::Origin, body2::Body)
    Ac = constraintmat(joint)
    An = nullspacemat(joint)
    Acᵀ = zerodimstaticadjoint(Ac)
    Anᵀ = zerodimstaticadjoint(An)
    Fτ = springforce(joint, body2.state) + damperforce(joint, body2.state)
    return Acᵀ * Ac * gc(joint, body2.state) + Anᵀ * An * Fτ
end

### Constraints and derivatives
## Discrete-time position wrappers (for dynamics)
gc(joint::ForceTranslational12, statea::State, stateb::State, Δt) = gc(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...) # USELESS
gc(joint::ForceTranslational12, stateb::State, Δt) = gc(joint, posargsnext(stateb, Δt)...) # USELESS
gc(joint::ForceTranslational12, statea::State, stateb::State) = gc(joint, posargsc(statea)..., posargsc(stateb)...)
gc(joint::ForceTranslational12, stateb::State) = gc(joint, posargsc(stateb)...)

@inline gc(joint::ForceTranslational12{T,N}) where {T,N} = szeros(T, N)


### Constraints and derivatives
## Position level constraints (for dynamics)
@inline function gc(joint::ForceTranslational12, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end
@inline function gc(joint::ForceTranslational12, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end


### Spring and damper
## Forces for dynamics
## Discrete-time position wrappers (for dynamics)
springforce(joint::ForceTranslational12, statea::State, stateb::State, Δt) = springforce(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...) #USELESS
springforce(joint::ForceTranslational12, stateb::State, Δt) = springforce(joint, posargsnext(stateb, Δt)...) #USELESS
springforce(joint::ForceTranslational12, statea::State, stateb::State) = springforce(joint, posargsc(statea)..., posargsc(stateb)...)
springforce(joint::ForceTranslational12, stateb::State) = springforce(joint, posargsc(stateb)...)

damperforce(joint::ForceTranslational12, statea::State, stateb::State, Δt) = damperforce(joint, statea.vsol[2], stateb.vsol[2]) #USELESS
damperforce(joint::ForceTranslational12, stateb::State, Δt) = damperforce(joint, stateb.vsol[2]) #USELESS
damperforce(joint::ForceTranslational12, statea::State, stateb::State) = damperforce(joint, statea.vsol[2], stateb.vsol[2])
damperforce(joint::ForceTranslational12, stateb::State) = damperforce(joint, stateb.vsol[2])

### Spring and damper
## Forces for dynamics
# Force applied by body a on body b expressed in world frame
@inline function springforce(joint::ForceTranslational12, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function springforce(joint::ForceTranslational12, xb::AbstractVector, qb::UnitQuaternion)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return force
end
# Force applied by body a on body b expressed in world frame
@inline function damperforce(joint::ForceTranslational12, va::AbstractVector, vb::AbstractVector)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * (vb - va)
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function damperforce(joint::ForceTranslational12, vb::AbstractVector)
    A = nullspacemat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vb
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return force
end


## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::ForceTranslational12, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    point2 = xb + vrotate(joint.vertices[2], qb)

    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * (Lmat(UnitQuaternion(point2)) - Lmat(UnitQuaternion(xa)))

    return X, Q
end
@inline function ∂g∂posb(joint::ForceTranslational12, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end
@inline function ∂g∂posb(joint::ForceTranslational12, xb::AbstractVector, qb::UnitQuaternion)
    X = I
    Q = 2 * VRᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))

    return X, Q
end

## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
@inline function ∂2g∂posaa(joint::ForceTranslational12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
    Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    XX = szeros(T, 9, 3)
    XQ = -kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) - kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = -kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Lpos)*∂L∂qsplit(T) + kron(VLᵀmat(qa)*Ltpos,2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::ForceTranslational12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = kron(VLᵀmat(qa),2*VLᵀmat(qa))*∂L∂qsplit(T)[:,SA[2; 3; 4]]
    QQ = kron(VLᵀmat(qa),2*VLᵀmat(qa)*Lmat(qb)*Lmat(UnitQuaternion(joint.vertices[2])))*∂Lᵀ∂qsplit(T) + kron(VLᵀmat(qa)*Lmat(qb)*Lᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa))*∂L∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::ForceTranslational12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = kron(Vmat(T),VLᵀmat(qa))*∂R∂qsplit(T) + kron(VRᵀmat(qa),Vmat(T))*∂Lᵀ∂qsplit(T)
    QX = szeros(T, 9, 3)
    QQ = kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb),2*VLᵀmat(qa))*∂R∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2]))*Rmat(qb)*Rᵀmat(qa),2*Vmat(T))*∂Lᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::ForceTranslational12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VLᵀmat(qa)*Rmat(qa)*Rᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*VLᵀmat(qa)*Rmat(qa))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::ForceTranslational12{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    QQ = kron(Vmat(T),2*VRᵀmat(qb)*Rmat(UnitQuaternion(joint.vertices[2])))*∂L∂qsplit(T) + kron(VLᵀmat(qb)*Rᵀmat(UnitQuaternion(joint.vertices[2])),2*Vmat(T))*∂Rᵀ∂qsplit(T)

    return XX, XQ, QX, QQ
end


### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::ForceTranslational12{T}, statea::State, stateb::State, clear::Bool) where T
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
@inline function applyFτ!(joint::ForceTranslational12{T}, stateb::State, clear::Bool) where T
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
@inline function ∂Fτ∂ua(joint::ForceTranslational12, statea::State, stateb::State)
    vertices = joint.vertices
    _, qa = posargsk(statea)

    BFa = -VLmat(qa) * RᵀVᵀmat(qa)
    Bτa = -skew(vertices[1])

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::ForceTranslational12, statea::State, stateb::State)
    vertices = joint.vertices
    xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)
    qbinvqa = qb\qa

    BFb = VLmat(qa) * RᵀVᵀmat(qa)
    Bτb = skew(vertices[2]) * VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::ForceTranslational12, stateb::State)
    vertices = joint.vertices
    _, qb = posargsk(stateb)

    BFb = I
    Bτb = skew(vertices[2]) * VLᵀmat(qb) * RVᵀmat(qb)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::ForceTranslational12{T}, statea::State, stateb::State) where T
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXa = szeros(T,3,3)
    FaQa = -2*VRᵀmat(qa)*Rmat(UnitQuaternion(F))*LVᵀmat(qa)
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,3)
    FbXa = szeros(T,3,3)
    FbQa = 2*VRᵀmat(qa)*Rmat(UnitQuaternion(F))*LVᵀmat(qa)
    τbXa = szeros(T,3,3)
    τbQa = 2*skew(vertices[2])*VLᵀmat(qb)*Rmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(F))*LVᵀmat(qa)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::ForceTranslational12{T}, statea::State, stateb::State) where T
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = 2*skew(vertices[2])*VLᵀmat(qb)*Lmat(qa)*Lmat(UnitQuaternion(F))*Lᵀmat(qa)*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::ForceTranslational12{T}, stateb::State) where T
    xb, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,3)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,3)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,3)
    τbXb = szeros(T,3,3)
    τbQb = 2*skew(vertices[2])*VLᵀmat(qb)*Lmat(UnitQuaternion(F))*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end


### Minimal coordinates
## Position and velocity offsets
@inline function getPositionDelta(joint::ForceTranslational12, body1::AbstractBody, body2::Body, x::SVector)
    Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
    return Δx
end
@inline function getVelocityDelta(joint::ForceTranslational12, body1::AbstractBody, body2::Body, v::SVector)
    Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
    return Δv
end

## Minimal coordinate calculation
@inline function minimalCoordinates(joint::ForceTranslational12, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
end
@inline function minimalCoordinates(joint::ForceTranslational12, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
end
@inline function minimalVelocities(joint::ForceTranslational12, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return nullspacemat(joint) * (stateb.vc - statea.vc)
end
@inline function minimalVelocities(joint::ForceTranslational12, body1::Origin, body2::Body)
    stateb = body2.state
    return nullspacemat(joint) * stateb.vc
end
