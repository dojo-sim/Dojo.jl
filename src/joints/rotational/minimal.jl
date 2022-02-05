@inline function get_position_delta(joint::Rotational, body1::Node, body2::Node,
		θ::SVector{N,T}) where {T,N}
    # axis angle representation
    θ = zerodimstaticadjoint(nullspace_mask(joint)) * θ
    # nθ = norm(θ)
    # if nθ == 0
    #     q = one(UnitQuaternion{T})
    # else
    #     q = UnitQuaternion(cos(nθ/2),(θ/nθ*sin(nθ/2))..., false)
    # end
	#
    # Δq = q * joint.qoffset # in body1 frame
    # return Δq
	Δq = axis_angle_to_quaternion(θ)
end

@inline function get_velocity_delta(joint::Rotational, body1::Node, body2::Node, ϕ::SVector)
    ϕ = zerodimstaticadjoint(nullspace_mask(joint)) * ϕ
    Δϕ = ϕ # in body1 frame
    return Δϕ
end

@inline function orientation_error(joint::Rotational, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    # (norm(Vmat(joint.qoffset)) > 1e-6) && (@warn "check the validity of this expression")

    q = inv(joint.qoffset) * inv(qa) * qb # rotation from frame b to frame a then joint_qoffset then spring_qoffset
    "q = Lᵀmat(qoffset) * Rmat(qb) * Tmat() * vector(qa)"
    "q = Lᵀmat(qoffset) * Lᵀmat(qa) * vector(qb)"
    # b --qb-> world --inv(qa)-> a --inv(qoffset)-> jointoffsetframe --inv(qoffset)-> springoffsetframe
    return q
end

function orientation_error_jacobian_configuration(jacobian_relative::Symbol,
        joint::Rotational, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion)
    (jacobian_relative == :parent) && (return Lᵀmat(joint.qoffset) * Rmat(qb) * Tmat() * LVᵀmat(qa))
    (jacobian_relative == :child) && (return Lᵀmat(joint.qoffset) * Lᵀmat(qa) * LVᵀmat(qb))
    return
end

################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::Rotational, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    q = inv(joint.qoffset) * inv(qa) * qb
    return nullspace_mask(joint) * rotation_vector(q)
end

@inline function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol,
        joint::Rotational{T}, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector,
        qb::UnitQuaternion) where T
    A = nullspace_mask(joint)
    q = inv(joint.qoffset) * inv(qa) * qb
    X = szeros(T,3,3)
    if jacobian_relative == :parent
        return A * [X ∂qrotation_vector(q) * Lᵀmat(joint.qoffset) * Rmat(qb) * Tmat() * LVᵀmat(qa)]
    elseif jacobian_relative == :child
        return A * [X ∂qrotation_vector(q) * Lᵀmat(joint.qoffset) * Lᵀmat(qa) * LVᵀmat(qb)]
    end
end

################################################################################
# Velocities
################################################################################
@inline function minimal_velocities(joint::Rotational, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15,
		stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
end

@inline function minimal_velocities(joint::Rotational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
	qoffset = joint.qoffset
    Δϕ = rotation_matrix(inv(qoffset)) * (vrotate(ϕb, inv(qa) * qb) - ϕa) # in offset frame
    return nullspace_mask(joint) * Δϕ
    # return nullspace_mask(joint) * (vrotate(ϕb, qa \ qb) - ϕa) # in body1's frame
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
        joint::Rotational, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_configuration_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_configuration_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
        joint::Rotational, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_velocity_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_velocity_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

@inline function minimal_velocities_jacobian_configuration_parent(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    X = szeros(T,3,3)
    Q = rotation_matrix(inv(qoffset)) * ∂qrotation_matrix_inv(qa, vrotate(ϕb, qb))
    ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qa)]
    return ∇xq
end
@inline function minimal_velocities_jacobian_configuration_child(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    X = szeros(T,3,3)
    Q = rotation_matrix(inv(qoffset) * inv(qa)) * ∂qrotation_matrix(qb, ϕb)
    ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qb)]
    return ∇xq
end
@inline function minimal_velocities_jacobian_velocity_parent(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    V = szeros(T,3,3)
    Ω = - rotation_matrix(inv(qoffset))
    ∇vϕ = nullspace_mask(joint) * [V Ω]
    return ∇vϕ
end
@inline function minimal_velocities_jacobian_velocity_child(joint::Rotational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
    qoffset = joint.qoffset
    V = szeros(T,3,3)
    Ω = rotation_matrix(inv(qoffset) * inv(qa) * qb)
    ∇vϕ = nullspace_mask(joint) * [V Ω]
    return ∇vϕ
end
