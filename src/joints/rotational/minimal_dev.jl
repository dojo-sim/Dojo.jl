@inline function get_position_delta(joint::Rotational, body1::Node, body2::Node, θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspace_mask(joint)) * θ
    nθ = norm(θ)
    if nθ == 0
        q = one(UnitQuaternion{T})
    else
        q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    end

    Δq = q * joint.qoffset # in body1 frame
    return Δq
end

@inline function get_velocity_delta(joint::Rotational, body1::Node, body2::Node, ω::SVector)
    ω = zerodimstaticadjoint(nullspace_mask(joint)) * ω
    Δω = ω # in body1 frame
    return Δω
end

@inline function minimal_coordinates(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.q2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
    q = qa \ qb / joint.qoffset
    return nullspace_mask(joint) * rotation_vector(q)
end

@inline function minimal_velocities(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.q2[1], statea.ϕ15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimal_velocities(joint::Rotational, qa::UnitQuaternion,
        ϕa::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    return nullspace_mask(joint) * (vrotate(ϕb, qa \ qb) - ϕa) # in body1's frame
end

function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol, joint::Rotational, qa::UnitQuaternion, qb::UnitQuaternion)
	c = minimal_coordinates(joint, qa, qb)
	jacobian_relative==:parent && (return [szeros(length(c), 3) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, UnitQuaternion(q..., false), qb), vector(qa))]) # * LVᵀmat(qa)]) 
	jacobian_relative==:child && (return [szeros(length(c), 3) FiniteDiff.finite_difference_jacobian(q -> minimal_coordinates(joint, qa, UnitQuaternion(q..., false)), vector(qb))]) # * LVᵀmat(qb)]) 
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
	joint::Rotational, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	v = minimal_velocities(joint, qa, ωa, qb, ωb)
	(jacobian_relative == :parent) && (return [szeros(length(v), 3) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, UnitQuaternion(q..., false), ωa, qb, ωb), vector(qa))]) # * LVᵀmat(qa)])
	(jacobian_relative == :child) && (return [szeros(length(v), 3) FiniteDiff.finite_difference_jacobian(q -> minimal_velocities(joint, qa, ωa, UnitQuaternion(q..., false), ωb), vector(qb))]) # * LVᵀmat(qb)])
	return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
	joint::Rotational, qa::UnitQuaternion, ωa::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	v = minimal_velocities(joint, qa, ωa, qb, ωb)
	(jacobian_relative == :parent) && (return [szeros(length(v), 3) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, qa, ω, qb, ωb), ωa)])
	(jacobian_relative == :child) && (return [szeros(length(v), 3) FiniteDiff.finite_difference_jacobian(ω -> minimal_velocities(joint, qa, ωa, qb, ω), ωb)])
	return
end
