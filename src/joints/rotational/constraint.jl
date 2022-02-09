mutable struct Rotational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N,Nb½}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    qoffset::UnitQuaternion{T} # in body1's frame

    spring
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
    unlimited_constraint(joint, xa, qa, xb, qb, η)
end
@inline function unlimited_constraint(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T}
    return constraint_mask(joint) * Vmat(orientation_error(joint, xa, qa, xb, qb)) # maybe we need to use rotation_vector instead of Vmat
end
#
# @inline function constraint_jacobian_configuration(joint::Rotational{T,Nλ,0,N}, η) where {T,Nλ,N}
#     return Diagonal(+1.00e-10 * sones(T,N))
# end

@inline function constraint_jacobian(jacobian_relative::Symbol, joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    Q = Vmat() * orientation_error_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb, attjac=false)
    return constraint_mask(joint) * [X Q]
end

@inline function constraint_jacobian_parent(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    # Q = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T)
    # return constraint_mask(joint) * [X Q]
    constraint_jacobian(:parent, joint, xa, qa, xb, qb, η)
end

@inline function constraint_jacobian_child(joint::Rotational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = szeros(T, 3, 3)
    # Q = VRᵀmat(joint.qoffset) * Lᵀmat(qa)
    # return constraint_mask(joint) * [X Q]
    constraint_jacobian(:child, joint, xa, qa, xb, qb, η)
end



################################################################################
# Impulse Transform
################################################################################
function impulse_transform_parent(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T}
    X = szeros(T, 3, 3)
    # QT = VRᵀmat(joint.qoffset) * Rmat(qb) * Tmat(T) * LVᵀmat(qa)
    Q = VLᵀmat(qa) * Tmat(T) * Rᵀmat(qb) * RVᵀmat(joint.qoffset)

    # τ = T(A->B)_a
    # Fτa = [Fa, τa] = [0, T(B->A)_a]
    # Q = -SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    # Q = -rotation_matrix(joint.qoffset)
    return [X; Q]
end

function impulse_transform_child(joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T}
    X = szeros(T, 3, 3)
    # QT = VRᵀmat(joint.qoffset) * Lᵀmat(qa) * LVᵀmat(qb)
    Q = VLᵀmat(qb) * Lmat(qa) * RVᵀmat(joint.qoffset)

    # τ = T(A->B)_a
    # Fτb = [Fb, τb] = [0, T(A->B)_b]
    # Q = rotation_matrix(inv(qb) * qa)
    # Q = rotation_matrix(inv(qb) * qa * joint.qoffset)
    return [X; Q]
end

################################################################################
 # Derivatives
################################################################################
function impulse_transform_parent_jacobian_parent(joint::Rotational{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapa'*p)/∂(xa,qa)
    Z3 = szeros(T,3,3)
    ∇Qqa = ∂qVLᵀmat(Tmat(T) * Rᵀmat(qb) * RVᵀmat(joint.qoffset) * p) * LVᵀmat(qa)
    return [Z3 Z3;
            Z3 ∇Qqa]
    return szeros(T,6,6)
    # error()
    # return 1e-10*SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
end

function impulse_transform_parent_jacobian_child(joint::Rotational{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapa'*p)/∂(xb,qb)
    Z3 = szeros(T,3,3)
    ∇Qqb = VLᵀmat(qa) * Tmat(T) * ∂qRᵀmat(RVᵀmat(joint.qoffset) * p) * LVᵀmat(qb)
    return [Z3 Z3;
            Z3 ∇Qqb]
    return szeros(T,6,6)
    # error()
    # return 1e-10*SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
end

function impulse_transform_child_jacobian_parent(joint::Rotational{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapb'*p)/∂(xa,qa)
    Z3 = szeros(T,3,3)
    ∇Qqa = VLᵀmat(qb) * ∂qLmat(RVᵀmat(joint.qoffset) * p) * LVᵀmat(qa)
    return [Z3 Z3;
            Z3 ∇Qqa]
    # Z3 = szeros(T,3,3)
    # ∇Qqa = rotation_matrix(inv(qb))  * ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
    # return [Z3 Z3;
    #         Z3 ∇Qqa]
end

function impulse_transform_child_jacobian_child(joint::Rotational{T,Nλ,0},
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapb'*p)/∂(xb,qb)
    Z3 = szeros(T,3,3)
    ∇Qqb = ∂qVLᵀmat(Lmat(qa) * RVᵀmat(joint.qoffset) * p) * LVᵀmat(qb)
    return [Z3 Z3;
            Z3 ∇Qqb]
    # Z3 = szeros(T,3,3)
    # ∇Qqb = ∂qrotation_matrix_inv(qb, rotation_matrix(qa)*p) * LVᵀmat(qb)
    # return [Z3 Z3;
    #         Z3 ∇Qqb]
end
