mutable struct Rotational{T,N,N̄} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T
    spring_offset::SVector{N̄,T}
    Fτ::SVector{3,T}

    function Rotational{T,N,N̄}(body1::Body, body2::Body;
            axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}),
            spring = zero(T), damper = zero(T), spring_offset = szeros(T,N̄)
        ) where {T,N,N̄}

        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)

        new{T,N,N̄}(V3, V12, qoffset, spring, damper, spring_offset, Fτ), body1.id, body2.id
    end
end

Rotational0 = Rotational{T,0,3} where T
Rotational1 = Rotational{T,1,2} where T
Rotational2 = Rotational{T,2,1} where T
Rotational3 = Rotational{T,3,0} where T

springforcea(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
springforceb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
springforceb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T, childid) where T = szeros(T, 6)

damperforcea(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
damperforceb(joint::Rotational{T,3}, body1::Body, body2::Body, Δt::T, childid) where T = szeros(T, 6)
damperforceb(joint::Rotational{T,3}, body1::Origin, body2::Body, Δt::T, childid) where T = szeros(T, 6)

# function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, constraint::Rotational{T,N}) where {T,N}
#     summary(io, constraint)
#     println(io,"")
#     println(io, " V3:      "*string(constraint.V3))
#     println(io, " V12:     "*string(constraint.V12))
#     println(io, " qoffset: "*string(constraint.qoffset))
# end

### Constraints and derivatives
## Position level constraints (for dynamics)
@inline function g(joint::Rotational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qa \ qb / joint.qoffset)))
    return Vmat(qa \ qb / joint.qoffset)
end

@inline function g(joint::Rotational, xb::AbstractVector, qb::UnitQuaternion)
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qb / joint.qoffset)))
    return Vmat(qb / joint.qoffset)
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)

    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T}, xb::AbstractVector, qb::UnitQuaternion) where T
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset)

    return X, Q
end

### Forcing
## Application of joint forces (for dynamics)
@inline function applyFτ!(joint::Rotational{T}, statea::State, stateb::State, Δt::T, clear::Bool) where T
    τ = joint.Fτ
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)

    τa = vrotate(-τ, qa) # in world coordinates
    τb = -τa # in world coordinates

    τa = vrotate(τa,inv(qa)) # in local coordinates
    τb = vrotate(τb,inv(qb)) # in local coordinates

    statea.τ2[end] += τa
    stateb.τ2[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

# for joints with origin as parent, torque is computed in child frame
@inline function applyFτ!(joint::Rotational{T}, stateb::State, Δt::T, clear::Bool) where T
    τ = joint.Fτ
    _, qb = posargs2(stateb)

    τb = τ
    τa = vrotate(-τb, qb) # from b frame to world coordinates

    stateb.τ2[end] += τb
    clear && (joint.Fτ = szeros(T,3))
    return
end

## Forcing derivatives (for linearization)
# Control derivatives
@inline function ∂Fτ∂ua(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    BFa = (szeros(T, 3, 3))
    Bτa = -I

    return [BFa; Bτa]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
    qbinvqa = qb \ qa

    BFb = (szeros(T, 3, 3))
    Bτb = VLmat(qbinvqa) * RᵀVᵀmat(qbinvqa)

    return [BFb; Bτb]
end
@inline function ∂Fτ∂ub(joint::Rotational{T}, stateb::State, Δt::T) where T
    _, qb = posargs2(stateb)

    BFb = (szeros(T, 3, 3))
    Bτb = I(3) # VLᵀmat(qb) * RVᵀmat(qb)

    return [BFb; Bτb]
end

# Position derivatives
@inline function ∂Fτ∂posa(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
    τ = joint.Fτ

    FaXa = szeros(T,3,3)
    FaQa = szeros(T,3,4)
    τaXa = szeros(T,3,3)
    τaQa = szeros(T,3,4)
    FbXa = szeros(T,3,3)
    FbQa = szeros(T,3,4)
    τbXa = szeros(T,3,3)
    τbQa = 2*VLᵀmat(qb)*Rmat(qb)*Rᵀmat(qa)*Rmat(UnitQuaternion(τ))#*LVᵀmat(qa)

    return FaXa, FaQa, τaXa, τaQa, FbXa, FbQa, τbXa, τbQa
end
@inline function ∂Fτ∂posb(joint::Rotational{T}, statea::State, stateb::State, Δt::T) where T
    _, qa = posargs2(statea)
    _, qb = posargs2(stateb)
    τ = joint.Fτ

    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = szeros(T,3,3)
    τaQb = szeros(T,3,4)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,4)
    τbXb = szeros(T,3,3)
    τbQb = 2*VLᵀmat(qb)*Lmat(qa)*Lmat(UnitQuaternion(τ))*Lᵀmat(qa)#*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end
@inline function ∂Fτ∂posb(joint::Rotational{T}, stateb::State, Δt::T) where T
    _, qb = posargs2(stateb)
    τ = joint.Fτ
    FaXb = szeros(T,3,3)
    FaQb = szeros(T,3,4)
    τaXb = szeros(T,3,3)
    τaQb = ∂vrotate∂q(-τ, qb)
    FbXb = szeros(T,3,3)
    FbQb = szeros(T,3,4)
    τbXb = szeros(T,3,3)
    τbQb = szeros(T,3,4)#2*VLᵀmat(qb)*Lmat(UnitQuaternion(τ))#*LVᵀmat(qb)

    return FaXb, FaQb, τaXb, τaQb, FbXb, FbQb, τbXb, τbQb
end

### Minimal coordinates
## Position and velocity offsets
@inline function getPositionDelta(joint::Rotational, body1::Body, body2::Body, θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspacemat(joint)) * θ
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end

    Δq = q * joint.qoffset # in body1 frame
    return Δq
end
@inline function getVelocityDelta(joint::Rotational, body1::Body, body2::Body, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = ω # in body1 frame
    # Δω = vrotate(ω, inv(body2.state.q2[1])*body1.state.q2[1]) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational, body1::Origin, body2::Body, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = ω # in body1 frame
    # Δω = vrotate(ω, inv(body2.state.q2[1])) # in body2 frame
    return Δω
end

## Minimal coordinate calculation (This could be directly calculated from g, but the rotation requires some special treatment)
@inline function minimalCoordinates(joint::Rotational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return minimalCoordinates(joint, statea.q2[1], stateb.q2[1])
end
@inline function minimalCoordinates(joint::Rotational, body1::Origin, body2::Body)
    stateb = body2.state
    return minimalCoordinates(joint, stateb.q2[1])
end
# useful for minimal to maximal coordinate mapping
@inline function minimalCoordinates(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    q = qa \ qb / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational, qb::UnitQuaternion)
    q = qb / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end
@inline function minimalCoordinates(joint::Rotational{T,0}, qb::UnitQuaternion) where {T}
    # @warn "removed 1st link quat"
    # q = qb / joint.qoffset
    # return vector(q)
    q = qb / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end

@inline function minimalVelocities(joint::Rotational, body1::Body, body2::Body)
    statea = body1.state
    stateb = body2.state
    return minimalVelocities(joint, statea.q2[1], statea.ϕ15, stateb.q2[1], stateb.ϕ15)
end
@inline function minimalVelocities(joint::Rotational, body1::Origin, body2::Body)
    stateb = body2.state
    return minimalVelocities(joint, stateb.q2[1], stateb.ϕ15)
end
@inline function minimalVelocities(joint::Rotational, qa::UnitQuaternion,
        ϕa::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    return nullspacemat(joint) * (vrotate(ϕb, qa \ qb) - ϕa) # in body1's frame
end
@inline function minimalVelocities(joint::Rotational, qb::UnitQuaternion, ϕb::AbstractVector)
    return nullspacemat(joint) * vrotate(ϕb, qb) # in body1's frame
end
