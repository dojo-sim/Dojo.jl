################################################################################
# GET: Coordinates Joint
################################################################################
@inline function minimal_coordinates(joint::Joint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_coordinates(joint, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
end

@inline function minimal_coordinates(joint::JointConstraint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
	Δx = minimal_coordinates(joint.translational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
    Δθ = minimal_coordinates(joint.rotational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
	return [Δx; Δθ]
end


################################################################################
# GET: Coordinates Joint
################################################################################
@inline function minimal_velocities_new(joint::JointConstraint, pnode::Node, cnode::Node, timestep)
	rot = joint.rotational
	tra = joint.translational
	pa = tra.vertices[1]
	pb = tra.vertices[2]
	qoffset = rot.qoffset
	Arot = nullspace_mask(rot)
	Atra = nullspace_mask(tra)

	xa, va, qa, ϕa = initial_configuration_velocity(pnode.state)
	xb, vb, qb, ϕb = initial_configuration_velocity(cnode.state)

	# Coordinates
	q = inv(qoffset) * inv(qa) * qb
	Δθ = Arot * rotation_vector(q)
	Δx = Atra * displacement(tra, xa, qa, xb, qb)

	# 1 step backward in time
	xa10 = next_position(xa, -va, timestep)
	qa10 = next_orientation(qa, -ϕa, timestep)
	xb10 = next_position(xb, -vb, timestep)
	qb10 = next_orientation(qb, -ϕb, timestep)

	# Previous step coordinates
	q10 = inv(qoffset) * inv(qa10) * qb10
	Δθ10 = Arot * rotation_vector(q10)
	Δx10 = Atra * displacement(tra, xa10, qa10, xb10, qb10)

	# Finite difference
	Δv = (Δx - Δx10) / timestep
    Δϕ = (Δθ - Δθ10) / timestep
	return [Δv; Δϕ]
end



################################################################################
# GET: Velocities Joint
################################################################################
@inline function minimal_velocities(joint::Joint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
    return minimal_velocities(joint, statea.x2[1], statea.v15, statea.q2[1], statea.ϕ15,
		stateb.x2[1], stateb.v15, stateb.q2[1], stateb.ϕ15)
end

function minimal_velocities_jacobian_configuration(jacobian_relative::Symbol,
        joint::Joint, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_configuration_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_configuration_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

function minimal_velocities_jacobian_velocity(jacobian_relative::Symbol,
        joint::Joint, xa::AbstractVector, va::AbstractVector,
        qa::UnitQuaternion, ϕa::AbstractVector, xb::AbstractVector,
        vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
    (jacobian_relative == :parent) && (return minimal_velocities_jacobian_velocity_parent(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    (jacobian_relative == :child) && (return minimal_velocities_jacobian_velocity_child(joint, xa, va, qa, ϕa, xb, vb, qb, ϕb))
    return
end

# ################################################################################
# # GET: Velocities Rotational
# ################################################################################
# @inline function minimal_velocities(joint::Rotational, xa::AbstractVector,
#         va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector)
# 	qoffset = joint.qoffset
#     Δϕ = rotation_matrix(inv(qoffset)) * (vrotate(ϕb, inv(qa) * qb) - ϕa) # in offset frame
#     return nullspace_mask(joint) * Δϕ
# end
#
# @inline function minimal_velocities_jacobian_configuration_parent(joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
#     qoffset = joint.qoffset
#     X = szeros(T,3,3)
#     Q = rotation_matrix(inv(qoffset)) * ∂qrotation_matrix_inv(qa, vrotate(ϕb, qb))
#     ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qa)]
#     return ∇xq
# end
# @inline function minimal_velocities_jacobian_configuration_child(joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
#     qoffset = joint.qoffset
#     X = szeros(T,3,3)
#     Q = rotation_matrix(inv(qoffset) * inv(qa)) * ∂qrotation_matrix(qb, ϕb)
#     ∇xq = nullspace_mask(joint) * [X Q*LVᵀmat(qb)]
#     return ∇xq
# end
# @inline function minimal_velocities_jacobian_velocity_parent(joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
#     qoffset = joint.qoffset
#     V = szeros(T,3,3)
#     Ω = - rotation_matrix(inv(qoffset))
#     ∇vϕ = nullspace_mask(joint) * [V Ω]
#     return ∇vϕ
# end
# @inline function minimal_velocities_jacobian_velocity_child(joint::Rotational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector) where T
#     qoffset = joint.qoffset
#     V = szeros(T,3,3)
#     Ω = rotation_matrix(inv(qoffset) * inv(qa) * qb)
#     ∇vϕ = nullspace_mask(joint) * [V Ω]
#     return ∇vϕ
# end


# ################################################################################
# # GET: Velocities Translational
# ################################################################################
# @inline function minimal_velocities(joint::Translational, xa::AbstractVector,
#         va::AbstractVector,  qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
# 	vertices = joint.vertices
#     pbcb_w = vrotate(-vertices[2], qb)
#     pbca_w = xa - (xb + vrotate(vertices[2], qb))
#     # Δvw = V(pb,B/A)w - V(pa,A/A)w
#     Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
#     Δv = vrotate(Δvw, inv(qa)) # in the a frame
#     return nullspace_mask(joint) * Δv
# end
#
# @inline function minimal_velocities_jacobian_configuration_parent(joint::Translational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
#     vertices = joint.vertices
#     pbcb_w = vrotate(-vertices[2], qb)
#     pbca_w = xa - (xb + vrotate(vertices[2], qb))
#     Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
#
#     X = - ∂pskew(vrotate(ωa, qa))
#     Q = - skew(pbca_w) * ∂qrotation_matrix(qa, ωa)
#     ∇xq = [X Q*LVᵀmat(qa)]
#     ∇xq = rotation_matrix(inv(qa)) * ∇xq + [szeros(T,3,3) ∂qrotation_matrix_inv(qa, Δvw) * LVᵀmat(qa)]
#     ∇xq = nullspace_mask(joint) * ∇xq
#     return ∇xq
# end
# @inline function minimal_velocities_jacobian_configuration_child(joint::Translational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
#     vertices = joint.vertices
#     pbcb_w = vrotate(-vertices[2], qb)
#
#     X = ∂pskew(vrotate(ωa, qa))
#     Q = ∂pskew(vrotate(ωb, qb)) * ∂qrotation_matrix(qb, -vertices[2]) + skew(pbcb_w) * ∂qrotation_matrix(qb, ωb)
#     Q += ∂pskew(vrotate(ωa, qa)) * ∂qrotation_matrix(qb, vertices[2])
#     ∇xq = [X Q*LVᵀmat(qb)]
#     ∇xq = rotation_matrix(inv(qa)) * ∇xq
#     ∇xq = nullspace_mask(joint) * ∇xq
#     return ∇xq
# end
# @inline function minimal_velocities_jacobian_velocity_parent(joint::Translational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
#     vertices = joint.vertices
#     pbca_w = xa - (xb + vrotate(vertices[2], qb))
#
#     V = - SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
#     Ω = - skew(pbca_w) * rotation_matrix(qa)
#     ∇vϕ = [V Ω]
#     ∇vϕ = rotation_matrix(inv(qa)) * ∇vϕ
#     ∇vϕ = nullspace_mask(joint) * ∇vϕ
#     return ∇vϕ
# end
# @inline function minimal_velocities_jacobian_velocity_child(joint::Translational{T},
#         xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
#         xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
#     vertices = joint.vertices
#     pbcb_w = vrotate(-vertices[2], qb)
#
#     V = SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
#     Ω = skew(pbcb_w) * rotation_matrix(qb)
#     ∇vϕ = [V Ω]
#     ∇vϕ = rotation_matrix(inv(qa)) * ∇vϕ
#     ∇vϕ = nullspace_mask(joint) * ∇vϕ
#     return ∇vϕ
# end
#
