mutable struct Rotational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    spring_type::Symbol # the rotational springs can be :sinusoidal or :linear, if linear then we need joint_limits to avoid the 180° singularity.
    Fτ::SVector{3,T}
end

function Rotational{T,Nλ}(body1::Node, body2::Node;
        axis::AbstractVector = szeros(T,3), qoffset::UnitQuaternion = one(UnitQuaternion{T}),
        spring = zero(T), damper = zero(T), spring_offset = szeros(T,3-Nλ),
        joint_limits = [szeros(T,0), szeros(T,0)],
        spring_type::Symbol = :sinusoidal,
    ) where {T,Nλ}

    V1, V2, V3 = orthogonal_rows(axis)
    V12 = [V1;V2]

    Fτ = zeros(T,3)
    Nb½ = length(joint_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    Rotational{T,Nλ,Nb,N,Nb½,N̄λ}(V3, V12, qoffset, spring, damper, spring_offset, joint_limits, spring_type, Fτ), body1.id, body2.id
end

Rotational0{T} = Rotational{T,0} where T
Rotational1{T} = Rotational{T,1} where T
Rotational2{T} = Rotational{T,2} where T
Rotational3{T} = Rotational{T,3} where T

# Position level constraints (for dynamics)
@inline function g(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    return constraintmat(joint) *  Vmat(qa \ qb / joint.qoffset)
end

@inline function ∂g∂z(joint::Rotational{T,Nλ,0,N}, η) where {T,Nλ,N}
    return Diagonal(+1.00e-10 * sones(T,N))
end

@inline function ∂g∂a(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return constraintmat(joint) * [X Q]
end

@inline function ∂g∂b(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    return constraintmat(joint) * [X Q]
end

@inline function Ga(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    Q = Q * LVᵀmat(qa)
    return constraintmat(joint) * [X Q]
end

@inline function Gb(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    Q = Q * LVᵀmat(qb)
    return constraintmat(joint) * [X Q]
end

## w/ Limits
@inline function g(joint::Rotational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    e1 = Vmat(qa \ qb / joint.qoffset)
    e2 = minimal_coordinates(joint, qa, qb)
    s, γ = get_sγ(joint, η)
    return [
            s .* γ;
            s[1:Nb½] - (joint.joint_limits[2] .- e2);
            s[Nb½ .+ (1:Nb½)] - (e2 .- joint.joint_limits[1]);
            constraintmat(joint) * e1;
           ]
end

@inline function ∂g∂z(joint::Rotational{T,Nλ,Nb,N}, η) where {T,Nλ,Nb,N}
    s, γ = get_sγ(joint, η)

    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline function ∂g∂a(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, UnitQuaternion(q..., false), xb, qb, η), vector(qa))
    return [X Q]
end

@inline function ∂g∂b(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, N, 3)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, UnitQuaternion(q..., false), η), vector(qb))
    return [X Q]
end

@inline function Ga(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return [
            zeros(Nb, 6);
            -nullspacemat(joint) * [X Q * LVᵀmat(qa)];
            nullspacemat(joint) * [X Q * LVᵀmat(qa)];
            constraintmat(joint) * [X Q * LVᵀmat(qa)];
           ]
end

@inline function Gb(joint::Rotational{T,Nλ,Nb,N}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    return [
            zeros(Nb, 6);
            -nullspacemat(joint) * [X Q * LVᵀmat(qb)];
            nullspacemat(joint) * [X Q * LVᵀmat(qb)];
            constraintmat(joint) * [X Q * LVᵀmat(qb)];
           ]
end

