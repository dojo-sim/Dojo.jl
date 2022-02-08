@inline function position_error(joint::Translational, xa::AbstractVector,
		qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)
	# TODO remove rotate
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
    rotate && (d = vrotate(d, inv(qa))) # in the a frame
    return d
end

function position_error_jacobian_configuration(jacobian_relative::Symbol,
        joint::Translational, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion)
    (jacobian_relative == :parent) && (return position_error_jacobian_configuration_parent(joint, xa, qa, xb, qb))
    (jacobian_relative == :child) && (return position_error_jacobian_configuration_child(joint, xa, qa, xb, qb))
    return
end

@inline function position_error_jacobian_configuration_parent(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
    X = - SMatrix{3,3,T,9}(Diagonal(sones(3)))
    Q = - ∂qrotation_matrix(qa, vertices[1])
    ∇xq = [X Q*LVᵀmat(qa)]
    ∇xq = rotation_matrix(inv(qa)) * ∇xq + [szeros(T,3,3) ∂qrotation_matrix_inv(qa, d) * LVᵀmat(qa)]
    return ∇xq
end

@inline function position_error_jacobian_configuration_child(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion) where T
    vertices = joint.vertices
    X = SMatrix{3,3,T,9}(Diagonal(sones(3)))
    Q = ∂qrotation_matrix(qb, vertices[2])
    ∇xq = [X Q*LVᵀmat(qb)]
    ∇xq = rotation_matrix(inv(qa)) * ∇xq
    return ∇xq
end

################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * position_error(joint, xa, qa, xb, qb)
end

@inline function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol, joint::Translational,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * position_error_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb)
end

################################################################################
# Velocities
################################################################################
@inline function minimal_velocities(joint::Translational, xa::AbstractVector,
        va::AbstractVector,  qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector)
	vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    # Δvw = V(pb,B/A)w - V(pa,A/A)w
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame
    Δv = vrotate(Δvw, inv(qa)) # in the a frame
    return nullspace_mask(joint) * Δv
end

@inline function minimal_velocities_jacobian_configuration_parent(joint::Translational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)
    pbca_w = xa - (xb + vrotate(vertices[2], qb))
    Δvw = vb + skew(pbcb_w) * vrotate(ωb, qb) - (va + skew(pbca_w) * vrotate(ωa, qa)) # in world frame

    X = - ∂pskew(vrotate(ωa, qa))
    Q = - skew(pbca_w) * ∂qrotation_matrix(qa, ωa)
    ∇xq = [X Q*LVᵀmat(qa)]
    ∇xq = rotation_matrix(inv(qa)) * ∇xq + [szeros(T,3,3) ∂qrotation_matrix_inv(qa, Δvw) * LVᵀmat(qa)]
    ∇xq = nullspace_mask(joint) * ∇xq
    return ∇xq
end
@inline function minimal_velocities_jacobian_configuration_child(joint::Translational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)

    X = ∂pskew(vrotate(ωa, qa))
    Q = ∂pskew(vrotate(ωb, qb)) * ∂qrotation_matrix(qb, -vertices[2]) + skew(pbcb_w) * ∂qrotation_matrix(qb, ωb)
    Q += ∂pskew(vrotate(ωa, qa)) * ∂qrotation_matrix(qb, vertices[2])
    ∇xq = [X Q*LVᵀmat(qb)]
    ∇xq = rotation_matrix(inv(qa)) * ∇xq
    ∇xq = nullspace_mask(joint) * ∇xq
    return ∇xq
end
@inline function minimal_velocities_jacobian_velocity_parent(joint::Translational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
    vertices = joint.vertices
    pbca_w = xa - (xb + vrotate(vertices[2], qb))

    V = - SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    Ω = - skew(pbca_w) * rotation_matrix(qa)
    ∇vϕ = [V Ω]
    ∇vϕ = rotation_matrix(inv(qa)) * ∇vϕ
    ∇vϕ = nullspace_mask(joint) * ∇vϕ
    return ∇vϕ
end
@inline function minimal_velocities_jacobian_velocity_child(joint::Translational{T},
        xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ωa::AbstractVector,
        xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ωb::AbstractVector) where T
    vertices = joint.vertices
    pbcb_w = vrotate(-vertices[2], qb)

    V = SMatrix{3,3,T,9}(Diagonal(sones(T,3)))
    Ω = skew(pbcb_w) * rotation_matrix(qb)
    ∇vϕ = [V Ω]
    ∇vϕ = rotation_matrix(inv(qa)) * ∇vϕ
    ∇vϕ = nullspace_mask(joint) * ∇vϕ
    return ∇vϕ
end
