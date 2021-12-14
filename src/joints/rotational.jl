mutable struct Rotational{T,N,N̄,Nl} <: Joint{T,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T
    spring_offset::SVector{N̄,T}
    joint_limits::Vector{SVector{Nl,T}} # lower and upper limits on the joint minimal coordinate angles
    Fτ::SVector{3,T}

    function Rotational{T,N,N̄}(body1::Component, body2::Component;
            axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}),
            spring = zero(T), damper = zero(T), spring_offset = szeros(T,N̄),
            joint_limits = [szeros(T,0), szeros(T,0)],
        ) where {T,N,N̄}

        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)
        nl = length(joint_limits[1])
        new{T,N,N̄,nl}(V3, V12, qoffset, spring, damper, spring_offset, joint_limits, Fτ), body1.id, body2.id
    end
end

joint_limits_length(joint::Rotational{T,N,N̄,Nl}) where {T,N,N̄,Nl} = Nl

Rotational0{T,Nl} = Rotational{T,0,3,Nl} where {T,Nl}
Rotational1{T,Nl} = Rotational{T,1,2,Nl} where {T,Nl}
Rotational2{T,Nl} = Rotational{T,2,1,Nl} where {T,Nl}
Rotational3{T,Nl} = Rotational{T,3,0,Nl} where {T,Nl}


# Position level constraints (for dynamics)
@inline function g(joint::Rotational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qa \ qb / joint.qoffset)))
    e = Vmat(qa \ qb / joint.qoffset)
    return constraintmat(joint) * e
end

@inline function g(joint::Rotational{T,N,N̄,0}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qb / joint.qoffset)))
    return constraintmat(joint) * Vmat(qb / joint.qoffset)
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return constraintmat(joint) * X, constraintmat(joint) * Q
end
@inline function ∂g∂posb(joint::Rotational{T,N,N̄,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)

    return constraintmat(joint) * X, constraintmat(joint) * Q
end
@inline function ∂g∂posb(joint::Rotational{T,N,N̄,0}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset)

    return constraintmat(joint) * X, constraintmat(joint) * Q
end


### w/ Limits

function get_sγ(joint::Rotational{T,N,N̄,Nl}, λ) where {T,N,N̄,Nl}
    s = λ[N .+ (1:(2 * Nl))] 
    γ = λ[N + 2 * Nl .+ (1:(2 * Nl))] 
    return s, γ
end

# Position level constraints (for dynamics)
@inline function g(joint::Rotational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    e1 = Vmat(qa \ qb / joint.qoffset)
    e2 = minimalCoordinates(joint, qa, qb)
    s, γ = get_sγ(joint, λ)
    return [
            constraintmat(joint) * e1;
            s[1:Nl] - (joint.joint_limits[2] - e2);
            s[Nl .+ (1:Nl)] - (e2 - joint.joint_limits[1]);
           ]
end

@inline function g(joint::Rotational{T,N,N̄,Nl}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    e1 = Vmat(qb / joint.qoffset)
    e2 = minimalCoordinates(joint, qb)
    s, γ = get_sγ(joint, λ)
    return [
            constraintmat(joint) * e1;
            s[1:Nl] - (joint.joint_limits[2] - e2);
            s[Nl .+ (1:Nl)] - (e2 - joint.joint_limits[1]);
           ]
end


## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = szeros(T, N + 2Nl + 2Nl, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, q, xb, qb, λ), qa) 
    return constraintmat(joint) * X, Q
end
@inline function ∂g∂posb(joint::Rotational{T,N,N̄,Nl}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = szeros(T, N + 2Nl + 2Nl, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, q, λ), qb) 
    return constraintmat(joint) * X, Q
end
@inline function ∂g∂posb(joint::Rotational{T,N,N̄,Nl}, xb::AbstractVector, qb::UnitQuaternion, λ) where {T,N,N̄,Nl}
    X = szeros(T, N + 2Nl + 2Nl, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xb, q, λ), qb)
    return constraintmat(joint) * X, Q
end

@inline function ∂g∂ʳself(joint::Rotational{T,N,N̄,0}, λ) where {T,N,N̄}
    return Diagonal(1e-10 * sones(T,N))
end

@inline function ∂g∂ʳself(joint::Rotational{T,N,N̄,Nl}, λ) where {T,N,N̄,Nl}
    # return 1e-10 * sones(T,N)
    [
     zeros(N, N + 4Nl);
     zeros(2Nl, N) Diagonal(ones(2Nl)) zeros(2Nl, 2Nl);
    ]
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
@inline function getPositionDelta(joint::Rotational, body1::Component, body2::Component, θ::SVector{N,T}) where {T,N}
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
@inline function getVelocityDelta(joint::Rotational, body1::Component, body2::Component, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = ω # in body1 frame
    # Δω = vrotate(ω, inv(body2.state.q2[1])*body1.state.q2[1]) # in body2 frame
    return Δω
end
@inline function getVelocityDelta(joint::Rotational, body1::Origin, body2::Component, ω::SVector)
    ω = zerodimstaticadjoint(nullspacemat(joint)) * ω
    Δω = ω # in body1 frame
    # Δω = vrotate(ω, inv(body2.state.q2[1])) # in body2 frame
    return Δω
end

## Minimal coordinate calculation (This could be directly calculated from g, but the rotation requires some special treatment)
@inline function minimalCoordinates(joint::Rotational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalCoordinates(joint, statea.q2[1], stateb.q2[1])
end
@inline function minimalCoordinates(joint::Rotational, body1::Origin, body2::Component)
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
    q = qb / joint.qoffset
    return nullspacemat(joint) * rotation_vector(q)
end

@inline function minimalVelocities(joint::Rotational, body1::Component, body2::Component)
    statea = body1.state
    stateb = body2.state
    return minimalVelocities(joint, statea.q2[1], statea.ϕ15, stateb.q2[1], stateb.ϕ15)
end
@inline function minimalVelocities(joint::Rotational, body1::Origin, body2::Component)
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
