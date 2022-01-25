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

@inline function constraint(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    return constraint_mask(joint) *  Vmat(qa \ qb / joint.qoffset)
end

@inline function constraint_jacobian_configuration(joint::Rotational{T,Nλ,0,N}, η) where {T,Nλ,N}
    return Diagonal(+1.00e-10 * sones(T,N))
end

@inline function constraint_jacobian_parent(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    return constraint_mask(joint) * [X Q]
end

@inline function constraint_jacobian_child(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    return constraint_mask(joint) * [X Q]
end

@inline function impulse_map_parent(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    Q = Q * LVᵀmat(qa)
    return constraint_mask(joint) * [X Q]
end

@inline function impulse_map_child(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    Q = Q * LVᵀmat(qb)
    return constraint_mask(joint) * [X Q]
end
