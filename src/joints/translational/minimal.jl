################################################################################
# Displacements
################################################################################
function displacement(joint::Translational, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion; rotate::Bool = true)

    vertices = joint.vertices
    d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa))
    rotate && (return vrotate(d, inv(qa))) : (return d)
end

function displacement_jacobian_configuration(relative::Symbol, joint::Translational{T}, 
    xa::AbstractVector, qa::UnitQuaternion, 
    xb::AbstractVector, qb::UnitQuaternion; attjac=true) where T

    vertices = joint.vertices

    if relative == :parent
        d = xb + vrotate(vertices[2], qb) - (xa + vrotate(vertices[1], qa)) # in the world frame
        X = -rotation_matrix(inv(qa))
        Q = -rotation_matrix(inv(qa)) * ∂qrotation_matrix(qa, vertices[1])
        Q += ∂qrotation_matrix_inv(qa, d)
        attjac && (Q *= LVᵀmat(qa))
    elseif relative == :child
        X = rotation_matrix(inv(qa))
        Q = rotation_matrix(inv(qa)) * ∂qrotation_matrix(qb, vertices[2])
        attjac && (Q *= LVᵀmat(qb))
    end

    return X, Q
end

################################################################################
# Coordinates
################################################################################
function minimal_coordinates(joint::Translational, xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion)
    return nullspace_mask(joint) * displacement(joint, xa, qa, xb, qb)
end

function minimal_coordinates_jacobian_configuration(relative::Symbol, joint::Translational,
        xa::AbstractVector, qa::UnitQuaternion, xb::AbstractVector, qb::UnitQuaternion; attjac::Bool=true)
    X, Q = displacement_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=attjac)
	return nullspace_mask(joint) * [X Q]
end

################################################################################
# Coordinates
################################################################################
function set_minimal_coordinates!(joint::Translational, 
    pnode::Node, cnode::Node,  
    timestep;
    Δx::AbstractVector=szeros(control_dimension(joint)))

    pa = joint.vertices[1]
    pb = joint.vertices[2]

    qa = pnode.state.q2[1]
    xa = pnode.state.x2[1]

    qb = cnode.state.q2[1]

    Aᵀ = zerodimstaticadjoint(nullspace_mask(joint))
    xb = xa + vrotate(pa + Aᵀ * Δx, qa) - vrotate(pb, qb)
    set_maximal_coordinates!(cnode; x=xb, q=cnode.state.q2[1])
    set_previous_configuration!(cnode, timestep)
    return nothing
end

################################################################################
# Velocities
################################################################################
function minimal_velocities(joint::Translational,
		xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
		xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
		timestep)
	A = nullspace_mask(joint)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)
	xb1 = next_position(xb, -vb, timestep)
	qb1 = next_orientation(qb, -ϕb, timestep)

	# Coordinates
	Δx = A * displacement(joint, xa, qa, xb, qb)
	# Previous step coordinates
	Δx1 = A * displacement(joint, xa1, qa1, xb1, qb1)

	# Finite difference
	Δv = (Δx - Δx1) / timestep
	return Δv
end

function minimal_velocities_jacobian_configuration(relative::Symbol, joint::Translational,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
    timestep)
    A = nullspace_mask(joint)

    # 1 step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)
    xb1 = next_position(xb, -vb, timestep)
    qb1 = next_orientation(qb, -ϕb, timestep)

    # Coordinates
    Δx = A * displacement(joint, xa, qa, xb, qb)
    # Previous step coordinates
    Δx1 = A * displacement(joint, xa1, qa1, xb1, qb1)

    # Finite difference
    # Δv = (Δx - Δx1) / timestep

    if relative == :parent
        X, Q = displacement_jacobian_configuration(:parent, joint, xa, qa, xb, qb, attjac=false)
        X1, Q1 = displacement_jacobian_configuration(:parent, joint, xa1, qa1, xb1, qb1, attjac=false)
        X1 *= -1.0
        Q1 *= -1.0 * Rmat(quaternion_map(-ϕa, timestep)) * timestep / 2
        Q *= LVᵀmat(qa) 
        Q1 *= LVᵀmat(qa)
        J = 1.0 / timestep * A * [X Q]
        J += 1.0 / timestep * A * [X1 Q1]
    elseif relative == :child 
        1.0 / timestep * (Δx - Δx1)
        X, Q = displacement_jacobian_configuration(:child, joint, xa, qa, xb, qb, attjac=false)
        X1, Q1 = displacement_jacobian_configuration(:child, joint, xa1, qa1, xb1, qb1, attjac=false)
        X1 *= -1.0
        Q1 *= -1.0 * Rmat(quaternion_map(-ϕb, timestep)) * timestep / 2
        Q *= LVᵀmat(qb) 
        Q1 *= LVᵀmat(qb)
        J = 1.0 / timestep * A * [X Q]
        J += 1.0 / timestep * A * [X1 Q1]
    end

    return J
end

function minimal_velocities_jacobian_velocity(relative::Symbol, joint::Translational,
    xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
    xb::AbstractVector, vb::AbstractVector, qb::UnitQuaternion, ϕb::AbstractVector,
    timestep)
    A = nullspace_mask(joint)

    # 1 step backward in time
    xa1 = next_position(xa, -va, timestep)
    qa1 = next_orientation(qa, -ϕa, timestep)
    xb1 = next_position(xb, -vb, timestep)
    qb1 = next_orientation(qb, -ϕb, timestep)

    # Coordinates
    Δx = A * displacement(joint, xa, qa, xb, qb)
    
    # Previous step coordinates
    Δx1 = A * displacement(joint, xa1, qa1, xb1, qb1)

    # Finite difference
    Δv = (Δx - Δx1) / timestep

    if relative == :parent
        X1, Q1 = displacement_jacobian_configuration(:parent, joint, xa1, qa1, xb1, qb1, attjac=false)
        X1 *= -timestep 
        Q1 *= -Lmat(qa) * quaternion_map_jacobian(-ϕa, timestep) * timestep / 2
        J = -1.0 / timestep * A * [X1 Q1]
    elseif relative == :child 
        X1, Q1 = displacement_jacobian_configuration(:child, joint, xa1, qa1, xb1, qb1, attjac=false)
        X1 *= -timestep 
        Q1 *= -Lmat(qb) * quaternion_map_jacobian(-ϕb, timestep) * timestep / 2
        J = -1.0 / timestep * A * [X1 Q1]
    end

    return J
end
