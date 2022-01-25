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

# Position level constraints (for dynamics)
@inline function constraint(joint::Translational{T,Nλ,0}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ}
    vertices = joint.vertices
    e = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    return constraint_mask(joint) * e
end

@inline function constraint_jacobian_configuration(joint::Translational{T,Nλ,0,N}, η) where {T,Nλ,N}
    return Diagonal(+1.00e-10 * sones(T,N))
end

## Derivatives NOT accounting for quaternion specialness
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

function impulse_map_parent(joint::Translational{T,Nλ,0}, statea::State, stateb::State, η, Δt) where {T,Nλ}
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

function impulse_map_child(joint::Translational{T,Nλ,0}, statea::State, stateb::State, η, Δt) where {T,Nλ}
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

### w/ Limits
# Position level constraints (for dynamics)
@inline function constraint(joint::Translational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    vertices = joint.vertices
    e1 = vrotate(xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)), inv(qa))
    e2 = minimal_coordinates(joint, xa, qa, xb, qb)

    s, γ = get_sγ(joint, η)

    return [
            s .* γ;
            s[1:Nb½] - (joint.joint_limits[2] - e2);
            s[Nb½ .+ (1:Nb½)] - (e2 - joint.joint_limits[1]);
            constraint_mask(joint) * e1;
           ]
end

@inline function constraint_jacobian_configuration(joint::Translational{T,Nλ,Nb,N}, η) where {T,Nλ,Nb,N}
    s, γ = get_sγ(joint, η)

    c1 = [Diagonal(γ + 1e-10 * sones(T, Nb)); Diagonal(sones(Nb)); szeros(Nλ, Nb)]
    c2 = [Diagonal(s + 1e-10 * sones(T, Nb)); szeros(Nb, Nb); szeros(Nλ, Nb)]
    c3 = [szeros(Nb, Nλ); szeros(Nb, Nλ); Diagonal(+1.00e-10 * sones(T, Nλ))]
    return [c1 c2 c3]
end

@inline function constraint_jacobian_parent(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    X = FiniteDiff.finite_difference_jacobian(x -> g(joint, x, qa, xb, qb, η), xa)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, UnitQuaternion(q..., false), xb, qb, η), vector(qa))
    return [X Q]
end

@inline function constraint_jacobian_child(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η)
    X = FiniteDiff.finite_difference_jacobian(x -> g(joint, xa, qa, x, qb, η), xb)
    Q = FiniteDiff.finite_difference_jacobian(q -> g(joint, xa, qa, xb, UnitQuaternion(q..., false), η), vector(qb))
    return [X Q]
end

function impulse_map_parent(joint::Translational, statea::State, stateb::State, η, Δt)
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_parent(joint, xa, qa, xb, qb, η)
end

function impulse_map_parent(joint::Translational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    X = -1.0 * transpose(rotation_matrix(qa))
    pb_a = rotation_matrix(inv(qa)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    ca_a = rotation_matrix(inv(qa)) * (xa) # body a com
    capb_a = pb_a - ca_a
    Q = - 1.0 * transpose(skew(capb_a))
    return [
            zeros(Nb, 6);
            -nullspace_mask(joint) * [X Q];
            nullspace_mask(joint) * [X Q];
            constraint_mask(joint) * [X Q];
           ]
end

function impulse_map_child(joint::Translational, statea::State, stateb::State, η, Δt)
    xa, qa = current_configuration(statea)
    xb, qb = current_configuration(stateb)
    impulse_map_child(joint, xa, qa, xb, qb, η)
end

function impulse_map_child(joint::Translational{T,Nλ,Nb,N,Nb½}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion, η) where {T,Nλ,Nb,N,Nb½}
    X = transpose(rotation_matrix(qa))
    pb_b = rotation_matrix(inv(qb)) * (xb + rotation_matrix(qb) * joint.vertices[2]) # body b kinematics point
    cb_b = rotation_matrix(inv(qb)) * (xb) # body b com
    cbpb_b = pb_b - cb_b
    Q = transpose(skew(cbpb_b) * rotation_matrix(inv(qb) * qa))
    return [
            zeros(Nb, 6);
            -nullspace_mask(joint) * [X Q];
            nullspace_mask(joint) * [X Q];
            constraint_mask(joint) * [X Q];
           ]
end

## Position and velocity offsets
@inline function get_position_delta(joint::Translational, body1::Node, body2::Node, x::SVector)
    Δx = zerodimstaticadjoint(nullspace_mask(joint)) * x # in body1 frame
    return Δx
end

@inline function get_velocity_delta(joint::Translational, body1::Node, body2::Node, v::SVector)
    Δv = zerodimstaticadjoint(nullspace_mask(joint)) * v # in body1 frame
    return Δv
end
