@inline function displacement(joint::Translational, xa::AbstractVector,
		qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    rotate && (return vrotate(d, inv(qa))) : (return d)
end

function displacement_jacobian_configuration(jacobian_relative::Symbol,
        joint::Translational, xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion; attjac=true)
    (jacobian_relative == :parent) && (return displacement_jacobian_configuration_parent(joint, xa, qa, xb, qb, attjac=attjac))
    (jacobian_relative == :child) && (return displacement_jacobian_configuration_child(joint, xa, qa, xb, qb, attjac=attjac))
    return
end

@inline function displacement_jacobian_configuration_parent(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac=true) where T
    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
    
    X = -rotation_matrix(inv(qa))
    Q = -rotation_matrix(inv(qa)) * ∂qrotation_matrix(qa, vertices[1])
    Q += ∂qrotation_matrix_inv(qa, d)
    attjac && (Q *= LVᵀmat(qa))
    
    return X, Q
end

@inline function displacement_jacobian_configuration_child(joint::Translational{T}, xa::AbstractVector,
        qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac=true) where T
    vertices = joint.vertices

    X = rotation_matrix(inv(qa))
    Q = rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, vertices[2])
    attjac && (Q *= LVᵀmat(qb))
    
    return X, Q
end

################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * displacement(joint, xa, qa, xb, qb)
end

@inline function minimal_coordinates_jacobian_configuration(jacobian_relative::Symbol, joint::Translational,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac::Bool=true)
    X, Q = displacement_jacobian_configuration(jacobian_relative, joint, xa, qa, xb, qb, attjac=attjac)
	return nullspace_mask(joint) * [X Q]
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
