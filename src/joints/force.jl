mutable struct Force12{T,N} <: FJoint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames

    spring::T
    damper::T

    Fτ::SVector{3,T}

    function Force12{T,N}(body1::AbstractBody, body2::AbstractBody;
            p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3), spring = zero(T), damper = zero(T)
        ) where {T,N}

        vertices = (p1, p2)
        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)

        new{T,N}(V3, V12, vertices, spring, damper, Fτ), body1.id, body2.id
    end
end

Force120 = Force12{T,0} where T
Force121 = Force12{T,1} where T
Force122 = Force12{T,2} where T
Force123 = Force12{T,3} where T

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Force12{T,N}) where {T,N}
    summary(io, constraint)
    println(io,"")
    println(io, " V3:       "*string(constraint.V3))
    println(io, " V12:      "*string(constraint.V12))
    println(io, " vertices: "*string(constraint.vertices))
end

### Constraints and derivatives
## Position level constraint wrappers
@inline function g(joint::Force12, body1::Body, body2::Body, Δt)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    Fτ = springforce(joint, body1.state, body2.state, Δt) + damperforce(joint, body1.state, body2.state)
    return Aᵀ * A * Fτ
end

@inline function g(joint::Force12, body1::Origin, body2::Body, Δt)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    Fτ = springforce(joint, body2.state, Δt) + damperforce(joint, body2.state)
    return Aᵀ * A * Fτ
end


## Position level constraints (for dynamics) in world frame
@inline function gc(joint::Force12, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
end
@inline function gc(joint::Force12, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return xb + vrotate(vertices[2], qb) - vertices[1]
end


### Spring and damper
## Forces for dynamics
## Discrete-time position wrappers (for dynamics)
springforce(joint::Force12, statea::State, stateb::State, Δt) = springforce(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)...)
springforce(joint::Force12, stateb::State, Δt) = springforce(joint, posargsnext(stateb, Δt)...)

damperforce(joint::Force12, statea::State, stateb::State) = damperforce(joint, statea.vsol[2], stateb.vsol[2])
damperforce(joint::Force12, stateb::State) = damperforce(joint, stateb.vsol[2])

### Spring and damper
## Forces for dynamics
# Force applied by body a on body b expressed in world frame
@inline function springforce(joint::Force12, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xa, qa, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function springforce(joint::Force12, xb::AbstractVector, qb::UnitQuaternion)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    distance = A * gc(joint, xb, qb)
    force = - Aᵀ * A * joint.spring * Aᵀ * distance  # Currently assumes same spring constant in all directions
    return force
end
# Force applied by body a on body b expressed in world frame
@inline function damperforce(joint::Force12, va::AbstractVector, vb::AbstractVector)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * (vb - va)
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return force
end
# Force applied by origin on body b expressed in world frame
@inline function damperforce(joint::Force12, vb::AbstractVector)
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    velocity = A * vb
    force = - Aᵀ * A * joint.damper * Aᵀ * velocity  # Currently assumes same damper constant in all directions
    return force
end

# Wrappers 2
∂g∂ʳposa(joint::Force12, statea::State, stateb::State, Δt) = ∂g∂ʳposa(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Force12, statea::State, stateb::State, Δt) = ∂g∂ʳposb(joint, posargsk(statea)..., posargsk(stateb)...)
∂g∂ʳposb(joint::Force12, stateb::State, Δt) = ∂g∂ʳposb(joint, posargsk(stateb)...)

# Derivatives accounting for quaternion specialness
# THIS IS USED TO INJECT THE SPRING-DAMPER FORCE INTO THE DYNAMICS, IT DOESN'T HAVE TO BE THE DERIVATIVE OF g WRT THE POS VARIABLES
@inline function ∂g∂ʳposa(joint::Force12{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    X = - Aᵀ * A # accounts for the fact that λsol[2] holds the force applied by body a on body b.
    Q = szeros(T, 3, 3)
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Force12{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    X = Aᵀ * A
    Q = szeros(T, 3, 3)
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Force12{T,N}, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    X = Aᵀ * A
    Q = szeros(T, 3, 3)
    return [X Q]
end

## Derivatives NOT accounting for quaternion specialness
# THIS IS USED IN DATAMAT, IT HAS TO BE THE DERIVATIVE OF g WRT THE POS VARIABLES
@inline function ∂g∂posa(joint::Force12{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    C = - Aᵀ * A * joint.spring * Aᵀ * A
    X = - C # xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    Q = - C * ∂vrotate∂q(joint.vertices[1], qa) # xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    return X, Q
end
@inline function ∂g∂posb(joint::Force12{T,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    C = - Aᵀ * A * joint.spring * Aᵀ * A
    X = C # xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    Q = C * ∂vrotate∂q(joint.vertices[2], qb) # xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    return X, Q
end
@inline function ∂g∂posb(joint::Force12{T,N}, xb::AbstractVector, qb::UnitQuaternion) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    C = - Aᵀ * A * joint.spring * Aᵀ * A
    X = C # xb + vrotate(vertices[2], qb)
    Q = C * ∂vrotate∂q(joint.vertices[2], qb) # xb + vrotate(vertices[2], qb)
    return X, Q
end

# Wrappers 2
∂g∂ʳvela(joint::Force12, statea::State, stateb::State, Δt) = ∂g∂ʳvela(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
∂g∂ʳvelb(joint::Force12, statea::State, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(statea)[2], posargsnext(statea, Δt)..., statea.vsol[2], statea.ωsol[2], posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)
∂g∂ʳvelb(joint::Force12, stateb::State, Δt) = ∂g∂ʳvelb(joint, posargsc(stateb)[2], posargsnext(stateb, Δt)..., stateb.vsol[2], stateb.ωsol[2], Δt)

# Derivatives accounting for quaternion specialness
@inline function ∂g∂ʳvela(joint::Force12{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
        vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    D = joint.damper * Aᵀ * A# damper
    S = Aᵀ * A * joint.spring * Aᵀ * A # spring
    V = D + S * Δt
    Ω = S * ∂vrotate∂q(joint.vertices[1], qa) * Lmat(q1a) * derivωbar(ωa, Δt) * Δt/2
    V *= Δt
    Ω *= Δt
    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Force12{T,N}, q1a::UnitQuaternion, xa::AbstractVector,
        qa::UnitQuaternion, va::AbstractVector, ωa::AbstractVector,
        q1b::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion,
        vb::AbstractVector, ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    D = - joint.damper * Aᵀ * A# damper
    S = - Aᵀ * A * joint.spring * Aᵀ * A # spring
    V = D + S * Δt
    Ω = S * ∂vrotate∂q(joint.vertices[2], qb) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    V *= Δt
    Ω *= Δt
    return [V Ω]
end
@inline function ∂g∂ʳvelb(joint::Force12{T,N}, q1b::UnitQuaternion, xb::AbstractVector,
        qb::UnitQuaternion, vb::AbstractVector,  ωb::AbstractVector, Δt) where {T,N}
    A = constraintmat(joint)
    Aᵀ = zerodimstaticadjoint(A)
    D = - joint.damper * Aᵀ * A# damper
    S = - Aᵀ * A * joint.spring * Aᵀ * A # spring
    V = D + S * Δt
    Ω = S * ∂vrotate∂q(joint.vertices[2], qb) * Lmat(q1b) * derivωbar(ωb, Δt) * Δt/2
    V *= Δt
    Ω *= Δt
    return [V Ω]
end

## vec(G) Jacobian (also NOT accounting for quaternion specialness in the second derivative: ∂(∂ʳg∂posx)∂y)
@inline function ∂2g∂posaa(joint::Force12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    Lpos = Lmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))
    Ltpos = Lᵀmat(UnitQuaternion(xb + vrotate(joint.vertices[2], qb) - xa))

    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)

    f = q -> ∂g∂ʳposa(joint, xa, UnitQuaternion(q...), xb, qb)[1:3, 4:6]
    df = ForwardDiff.jacobian(f, [qa.w; qa.x; qa.y; qa.z])
    # @show df

    QQ = df#szeros(T, 9, 4) 
    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posab(joint::Force12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)

    f = q -> ∂g∂ʳposa(joint, xa, qa, xb, UnitQuaternion(q...))[1:3, 4:6]
    df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])
    # @show df
    QQ = df#szeros(T, 9, 4) 

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posba(joint::Force12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)

    f = q -> ∂g∂ʳposb(joint, xa, UnitQuaternion(q...), xb, qb)[1:3, 4:6]
    df = ForwardDiff.jacobian(f, [qa.w; qa.x; qa.y; qa.z])
    # @show df

    QQ = df#szeros(T, 9, 4)
    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Force12{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)

    f = q -> ∂g∂ʳposb(joint, xa, qa, xb, UnitQuaternion(q...))[1:3, 4:6]
    df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])
    # @show df

    QQ = df#szeros(T, 9, 4)

    return XX, XQ, QX, QQ
end
@inline function ∂2g∂posbb(joint::Force12{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    XX = szeros(T, 9, 3)
    XQ = szeros(T, 9, 4)
    QX = szeros(T, 9, 3)
    
    f = q -> ∂g∂ʳposb(joint, xb, UnitQuaternion(q...))[1:3, 4:6]
    df = ForwardDiff.jacobian(f, [qb.w; qb.x; qb.y; qb.z])

    # @show df
    QQ = df#szeros(T, 9, 4)
    return XX, XQ, QX, QQ
end

### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Force12{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
    F = joint.Fτ
    vertices = joint.vertices
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)

    Fa = -F
    Fb =  F

    τa = vrotate(torqueFromForce(Fa, vrotate(vertices[1], qa)),inv(qa)) # in local coordinates R(qainv) * skew(R(qa) * vert) * Fa
    τb = vrotate(torqueFromForce(Fb, vrotate(vertices[2], qb)),inv(qb)) # in local coordinates

    statea.Fk[end] += Fa
    statea.τk[end] += τa
    stateb.Fk[end] += Fb
    stateb.τk[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end
@inline function applyFτ!(joint::Force12{T}, stateb::State, Δt::T, clear::Bool) where T
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
@inline function ∂Fτ∂ua(joint::Force12, statea::State, stateb::State)
    vertices = joint.vertices
    _, qa = posargsk(statea)

    BFa = -I(3)
    Bτa = -skew(vertices[1]) * VLmat(qa) * RᵀVᵀmat(qa)

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Force12, statea::State, stateb::State)
    vertices = joint.vertices
    # xa, qa = posargsk(statea)
    xb, qb = posargsk(stateb)

    BFb = I(3)
    Bτb = skew(vertices[2]) * VLmat(qb) * RᵀVᵀmat(qb)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Force12, stateb::State)
    vertices = joint.vertices
    _, qb = posargsk(stateb)

    BFb = I
    Bτb = skew(vertices[2]) * VLᵀmat(qb) * RVᵀmat(qb)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Force12{T}, statea::State, stateb::State) where T
    _, qa = posargsk(statea)
    _, qb = posargsk(stateb)
    F = joint.Fτ
    vertices = joint.vertices

    FaXa = szeros(T,3,3)
    FaQa = szeros(T,3,3)
    τaXa = szeros(T,3,3)
    τaQa = 2*skew(vertices[1])*VLᵀmat(qa)*Lmat(UnitQuaternion(F))*LVᵀmat(qa)
    FbXa = szeros(T,3,3)
    FbQa = szeros(T,3,3)
    τbXa = szeros(T,3,3)
    τbQa = szeros(T,3,3)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Force12{T}, statea::State, stateb::State) where T
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
    τbQb = 2*skew(vertices[2])*VLᵀmat(qb)*Lmat(UnitQuaternion(F))*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Force12{T}, stateb::State) where T
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


# ### Minimal coordinates
# ## Position and velocity offsets
# @inline function getPositionDelta(joint::Force12, body1::AbstractBody, body2::Body, x::SVector)
#     Δx = zerodimstaticadjoint(nullspacemat(joint)) * x # in body1 frame
#     return Δx
# end
# @inline function getVelocityDelta(joint::Force12, body1::AbstractBody, body2::Body, v::SVector)
#     Δv = zerodimstaticadjoint(nullspacemat(joint)) * v # in body1 frame
#     return Δv
# end

# ## Minimal coordinate calculation
# @inline function minimalCoordinates(joint::Force12, body1::Body, body2::Body)
#     statea = body1.state
#     stateb = body2.state
#     return nullspacemat(joint) * g(joint, statea.xc, statea.qc, stateb.xc, stateb.qc)
# end
# @inline function minimalCoordinates(joint::Force12, body1::Origin, body2::Body)
#     stateb = body2.state
#     return nullspacemat(joint) * g(joint, stateb.xc, stateb.qc)
# end
# @inline function minimalVelocities(joint::Force12, body1::Body, body2::Body)
#     statea = body1.state
#     stateb = body2.state
#     return nullspacemat(joint) * (stateb.vc - statea.vc)
# end
# @inline function minimalVelocities(joint::Force12, body1::Origin, body2::Body)
#     stateb = body2.state
#     return nullspacemat(joint) * stateb.vc
# end
# nullspacemat(force1)
