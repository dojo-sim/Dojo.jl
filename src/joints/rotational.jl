mutable struct Rotational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    Fτ::SVector{3,T}

    function Rotational{T,Nλ}(body1::Component, body2::Component;
            axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}),
            spring = zero(T), damper = zero(T), spring_offset = szeros(T,3-Nλ),
            joint_limits = [szeros(T,0), szeros(T,0)],
        ) where {T,Nλ}

        V1, V2, V3 = orthogonalrows(axis)
        V12 = [V1;V2]

        Fτ = zeros(T,3)
        Nb½ = length(joint_limits[1])
        Nb = 2Nb½
        N̄λ = 3 - Nλ
        N = Nλ + 2Nb
        new{T,Nλ,Nb,N,Nb½,N̄λ}(V3, V12, qoffset, spring, damper, spring_offset, joint_limits, Fτ), body1.id, body2.id
    end
end

Rotational0{T} = Rotational{T,0} where T
Rotational1{T} = Rotational{T,1} where T
Rotational2{T} = Rotational{T,2} where T
Rotational3{T} = Rotational{T,3} where T


# Position level constraints (for dynamics)
@inline function g(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qa \ qb / joint.qoffset)))
    return constraintmat(joint) *  Vmat(qa \ qb / joint.qoffset)
end

@inline function g(joint::Rotational{T,Nλ,0}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    # typeof(joint) <: Rotational{Float64,1} && println(scn.(Vmat(qb / joint.qoffset)))
    return constraintmat(joint) * Vmat(qb / joint.qoffset)
end

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return constraintmat(joint) * X, constraintmat(joint) * Q
end
@inline function ∂g∂posb(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    return constraintmat(joint) * X, constraintmat(joint) * Q
end
@inline function ∂g∂posb(joint::Rotational{T,Nλ,0}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset)
    return constraintmat(joint) * X, constraintmat(joint) * Q
end

@inline function ∂g∂ʳposa(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X, Q = ∂g∂posa(joint, xa, qa, xb, qb, η)
    Q = Q * LVᵀmat(qa)
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X, Q = ∂g∂posb(joint, xa, qa, xb, qb, η)
    Q = Q * LVᵀmat(qb)
    return [X Q]
end
@inline function ∂g∂ʳposb(joint::Rotational{T,Nλ,0}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X, Q = ∂g∂posb(joint, xb, qb, η)
    Q = Q * LVᵀmat(qb)
    return [X Q]
end

### w/ Limits
# Position level constraints (for dynamics)
@inline function g(joint::Rotational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = Vmat(qa \ qb / joint.qoffset)
    e2 = minimalCoordinates(joint, qa, qb)
    s, γ = get_sγ(joint, η)
    return [
            s[1:Nb½] - (joint.joint_limits[2] .- e2);
            s[Nb½ .+ (1:Nb½)] - (e2 .- joint.joint_limits[1]);
            constraintmat(joint) * e1;
           ]
end

@inline function g(joint::Rotational{T,Nλ,Nb,N,Nb½}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = Vmat(qb / joint.qoffset)
    e2 = minimalCoordinates(joint, qb)
    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[1:Nb½] - (joint.joint_limits[2] - e2);
            s[Nb½ .+ (1:Nb½)] - (e2 - joint.joint_limits[1]);
            constraintmat(joint) * e1;
           ]
end


## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂posa(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, UnitQuaternion(q..., false), xb, qb, η), vector(qa))
    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, UnitQuaternion(q..., false), η), vector(qb))
    return X, Q
end
@inline function ∂g∂posb(joint::Rotational{T,Nλ,Nb,N}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xb, UnitQuaternion(q..., false), η), vector(qb))
    return X, Q
end

@inline function ∂g∂ʳposa(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return [
            zeros(Nb, 6);
            -nullspacemat(joint) * [X Q * LVᵀmat(qa)];
             nullspacemat(joint) * [X Q * LVᵀmat(qa)];
            constraintmat(joint) * [X Q * LVᵀmat(qa)];
           ]
end
@inline function ∂g∂ʳposb(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    return [
            zeros(Nb, 6);
            -nullspacemat(joint) * [X Q * LVᵀmat(qb)];
             nullspacemat(joint) * [X Q * LVᵀmat(qb)];
            constraintmat(joint) * [X Q * LVᵀmat(qb)];
           ]
end
@inline function ∂g∂ʳposb(joint::Rotational{T,Nλ,Nb,N}, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset)
    return [
            zeros(Nb, 6);
            -nullspacemat(joint) * [X Q * LVᵀmat(qb)];
             nullspacemat(joint) * [X Q * LVᵀmat(qb)];
            constraintmat(joint) * [X Q * LVᵀmat(qb)];
           ]
end


@inline function ∂g∂ʳself(joint::Rotational{T,Nλ,0,N}, η) where {T,Nλ,N}
    return Diagonal(1e-10 * sones(T,N))
end

@inline function ∂g∂ʳself(joint::Rotational{T,Nλ,Nb,N}, η) where {T,Nλ,Nb,N}
    s, γ = get_sγ(joint, η)

    [
     Diagonal(γ) Diagonal(s) zeros(Nb, Nλ);
     Diagonal(ones(Nb)) zeros(Nb, Nb + Nλ);
     zeros(Nλ, N);
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
