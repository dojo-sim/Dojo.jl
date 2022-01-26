mutable struct Translational{T,Nλ,Nb,N,Nb½,N̄λ} <: Joint{T,Nλ,Nb,N}
    V3::Adjoint{T,SVector{3,T}} # in body1's frame
    V12::SMatrix{2,3,T,6} # in body1's frame
    vertices::NTuple{2,SVector{3,T}} # in body1's & body2's frames
    spring::T
    damper::T
    spring_offset::SVector{N̄λ,T}
    joint_limits::Vector{SVector{Nb½,T}} # lower and upper limits on the joint minimal coordinate angles
    spring_type::Symbol # the rotational springs can be :sinusoidal or :linear, if linear then we need joint_limits to avoid the 180° singularity.
    Fτ::SVector{3,T}
end

function Translational{T,Nλ}(body1::Node, body2::Node;
        p1::AbstractVector = szeros(T,3), p2::AbstractVector = szeros(T,3), axis::AbstractVector = szeros(T,3),
        spring = zero(T), damper = zero(T), spring_offset = szeros(T,3-Nλ),
        joint_limits = [szeros(T,0), szeros(T,0)],
        spring_type::Symbol = :sinusoidal,
    ) where {T,Nλ}
    vertices = (p1, p2)
    V1, V2, V3 = orthogonal_rows(axis)
    V12 = [V1;V2]
    Fτ = zeros(T,3)
    Nb½ = length(joint_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    Translational{T,Nλ,Nb,N,Nb½,N̄λ}(V3, V12, vertices, spring, damper, spring_offset, joint_limits, spring_type, Fτ), body1.id, body2.id
end

Translational0{T} = Translational{T,0} where T
Translational1{T} = Translational{T,1} where T
Translational2{T} = Translational{T,2} where T
Translational3{T} = Translational{T,3} where T

@inline function position_error(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    vertices = joint.vertices
    return vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
end

@inline function constraint(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    vertices = joint.vertices
    e = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    return constraint_mask(joint) * e
end

@inline function constraint_jacobian_configuration(joint::Translational{T,Nλ,0,N}, η) where {T,Nλ,N}
    return Diagonal(+1.00e-10 * sones(T,N))
end

@inline function constraint_jacobian_parent(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    point2 = xb + vrotate(joint.vertices[2], qb)
    X = -VLᵀmat(qa) * RVᵀmat(qa)
    Q = ∂vrotate∂q(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * Tmat()
    Q += ∂vrotate∂p(point2 - (xa + vrotate(joint.vertices[1], qa)), inv(qa)) * -∂vrotate∂q(joint.vertices[1], qa)
    return constraint_mask(joint) * [X Q]
end

@inline function constraint_jacobian_child(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = VLᵀmat(qa) * RVᵀmat(qa)
    Q = 2 * VLᵀmat(qa) * Rmat(qa) * Rᵀmat(qb) * Rmat(UnitQuaternion(joint.vertices[2]))
    return constraint_mask(joint) * [X Q]
end

function impulse_map_parent(joint::Translational{T,Nλ,0}, statea::State, stateb::State, η, timestep) where {T,Nλ}
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_parent(joint, xa, qa, xb, qb, η)
end

function impulse_map_parent(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = -1.0 * transpose(rotation_matrix(qa))
    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com
    capb_a = pb_a - ca_a
    Q = - 1.0 * transpose(skew(capb_a))
    return constraint_mask(joint) * [X Q]
end

function impulse_map_child(joint::Translational{T,Nλ,0}, statea::State, stateb::State, η, timestep) where {T,Nλ}
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_child(joint, xa, qa, xb, qb, η)
end

function impulse_map_child(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    X = transpose(rotation_matrix(qa))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb) * qa))
    return constraint_mask(joint) * [X Q]
end

@inline function get_position_delta(joint::Translational, body1::Node, body2::Node, x::SVector)
    Δx = zerodimstaticadjoint(nullspace_mask(joint)) * x # in body1 frame
    return Δx
end

@inline function get_velocity_delta(joint::Translational, body1::Node, body2::Node, v::SVector)
    Δv = zerodimstaticadjoint(nullspace_mask(joint)) * v # in body1 frame
    return Δv
end

################################################################################
# Force Mapping & Projector
################################################################################

function force_mapa(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T}
    X = -1.0 * rotation_matrix(qa)
    # pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    # ca_a = rotation_matrix(inv(qa)) * (xa) # body a com
    # capb_a = pb_a - ca_a
    # Q = - 1.0 * skew(capb_a)

    capb_a = rotation_matrix(inv(qa)) * (xb - xa + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    Q = - 1.0 * skew(capb_a)
    return [X; Q]
end

function force_mapb(joint::Translational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where {T}
    X = rotation_matrix(qa)
    # pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    # cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    # cbpb_b = pb_b - cb_b
    # Q = skew(cbpb_b) * rotation_matrix(inv(qb) * qa)

    cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
    Q = rotation_matrix(inv(qb)) * skew(cbpb_w) * rotation_matrix(qa)
    return [X; Q]
end


################################################################################
 # Derivatives
################################################################################

function ∂aforce_mapb(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapb'*p)/∂(xa,qa)
    Z3 = szeros(T,3,3)
    cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point

    ∇Xqa = ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
    ∇Qqa = rotation_matrix(inv(qb)) * skew(cbpb_w) * ∂qrotation_matrix(qa, p) * LVᵀmat(qa)
    return [Z3 ∇Xqa;
            Z3 ∇Qqa]
end

function ∂bforce_mapb(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapb'*p)/∂(xb,qb)
    Z3 = szeros(T,3,3)
    cbpb_w = rotation_matrix(qb) * joint.vertices[2] # body b kinematics point
    ∇Qqb = ∂qrotation_matrix_inv(qb, skew(cbpb_w) * rotation_matrix(qa) * p)
    ∇Qqb += rotation_matrix(inv(qb)) * ∂pskew(rotation_matrix(qa) * p) * ∂qrotation_matrix(qb, joint.vertices[2])
    ∇Qqb *= LVᵀmat(qb)
    return [Z3 Z3;
            Z3 ∇Qqb]
end

function ∂aforce_mapa(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapa'*p)/∂(xa,qa)
    Z3 = szeros(T,3,3)

    ∇Xqa = -∂qrotation_matrix(qa, p) * LVᵀmat(qa)
    ∇Qxa =  ∂pskew(p) * rotation_matrix(inv(qa))
    ∇Qqa = -∂pskew(p) * ∂qrotation_matrix_inv(qa, xb - xa + rotation_matrix(qb) * joint.vertices[2]) * LVᵀmat(qa)
    return [Z3   ∇Xqa;
            ∇Qxa ∇Qqa]
end

function ∂bforce_mapa(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, p) where {T,Nλ}
    # ∂(force_mapa'*p)/∂(xb,qb)
    Z3 = szeros(T,3,3)

    ∇Qxb = -∂pskew(p) * rotation_matrix(inv(qa))
    ∇Qqb = -∂pskew(p) * rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, joint.vertices[2]) * LVᵀmat(qb)
    return [Z3   Z3;
            ∇Qxb ∇Qqb]
end
