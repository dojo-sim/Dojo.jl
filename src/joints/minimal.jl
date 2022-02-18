################################################################################
# Coordinates
################################################################################
@inline function minimal_coordinates(joint::JointConstraint, body1::Node, body2::Node)
    statea = body1.state
    stateb = body2.state
	Δx = minimal_coordinates(joint.translational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
    Δθ = minimal_coordinates(joint.rotational, statea.x2[1], statea.q2[1], stateb.x2[1], stateb.q2[1])
	return [Δx; Δθ]
end

function set_minimal_coordinates!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
        Δx::AbstractVector=szeros(control_dimension(joint.translational)),
        Δθ::AbstractVector=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since xb = fct(qb, Δx)
    set_minimal_coordinates!(pnode, cnode, joint.rotational, timestep; Δθ=Δθ)
    set_minimal_coordinates!(pnode, cnode, joint.translational, timestep; Δx=Δx)
    return nothing
end

function set_position!(mechanism, joint::JointConstraint, xθ; iter::Bool=true)
    if !iter
        set_joint_position!(mechanism, joint, xθ)
    else
        currentvals = minimal_coordinates(mechanism)
        set_joint_position!(mechanism, joint, xθ)
        for id in recursivedirectchildren!(mechanism.system, joint.id)
            node = get_node(mechanism, id)
            if node isa JointConstraint
                set_joint_position!(mechanism, node, currentvals[id])
            end
        end
    end
    return
end

function set_joint_position!(mechanism, joint::JointConstraint{T,N,Nc}, xθ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += joint_length(element)
    end
    @assert length(xθ)==3*Nc-Nλ

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δx = xθ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δθ = xθ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_coordinates!(body1, body2, joint, mechanism.timestep, Δx=Δx, Δθ=Δθ)
    return body2.state.x2[1], body2.state.q2[1]
end

################################################################################
# Velocities
################################################################################
@inline function minimal_coordinates_velocities(joint::JointConstraint,
    pnode::Node, cnode::Node, timestep)
    Δxθ = minimal_coordinates(joint, pnode, cnode)
    Δvϕ = minimal_velocities(joint, pnode, cnode, timestep)
    return [Δxθ; Δvϕ]
end

@inline function minimal_velocities(joint::JointConstraint, pnode::Node, cnode::Node, timestep)
    Δv = minimal_velocities(joint.translational, pnode, cnode, timestep)
    Δϕ = minimal_velocities(joint.rotational, pnode, cnode, timestep)
    return [Δv; Δϕ]
end

function set_minimal_velocities!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
        Δv=szeros(control_dimension(joint.translational)),
        Δϕ=szeros(control_dimension(joint.rotational)))
    # We need to set the minimal coordinates of the rotational joint first
    # since vb = fct(ϕb, Δv)
    # set_minimal_velocities!(pnode, cnode, joint.rotational; Δϕ=Δϕ)
    # set_minimal_velocities!(pnode, cnode, joint.translational; Δv=Δv)
	vb, ϕb = set_minimal_velocities(joint, initial_configuration_velocity(pnode.state)...,
	 	current_configuration(cnode.state)..., timestep, Δv=Δv, Δϕ=Δϕ)
	set_velocity!(cnode; v=vb, ω=ϕb)
	set_previous_configuration!(cnode, timestep)
    return nothing
end

function set_minimal_velocities(joint::JointConstraint,
		xa::AbstractVector, va::AbstractVector, qa::UnitQuaternion, ϕa::AbstractVector,
		xb::AbstractVector, qb::UnitQuaternion, timestep;
        Δv=szeros(control_dimension(joint.translational)),
        Δϕ=szeros(control_dimension(joint.rotational)))
	rot = joint.rotational
	tra = joint.translational
	pa = tra.vertices[1]
    pb = tra.vertices[2]
	qoffset = rot.qoffset
	Arotᵀ = zerodimstaticadjoint(nullspace_mask(rot))
	Atraᵀ = zerodimstaticadjoint(nullspace_mask(tra))

	Δx = minimal_coordinates(joint.translational, xa, qa, xb, qb)

	# 1 step backward in time
	xa1 = next_position(xa, -va, timestep)
	qa1 = next_orientation(qa, -ϕa, timestep)

	# Finite difference
	Δx1 = Δx .- Δv * timestep

	# Δθ is expressed in along the joint's nullspace axes, in pnode's offset frame
	Δq2 = inv(qoffset) * qa * qb
	Δq1 = Δq2 * inv(axis_angle_to_quaternion(Arotᵀ * Δϕ * timestep))
	qb1 = qa1 * qoffset * Δq1

    xb1 = xa1 + vrotate(pa + Atraᵀ*Δx1, qa1) - vrotate(pb, qb1)

	# Finite difference
	vb = (xb - xb1) / timestep
	ϕb = angular_velocity(qb1, qb, timestep)
    return vb, ϕb
end

function set_velocity!(mechanism, joint::JointConstraint{T,N,Nc}, vϕ) where {T,N,Nc}
    Nλ = 0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        Nλ += joint_length(element)
    end

    # bodies
    body1 = get_body(mechanism, joint.parent_id)
    body2 = get_body(mechanism, joint.child_id)

    Δv = vϕ[SUnitRange(joint.minimal_index[1][1], joint.minimal_index[1][2])]
    Δϕ = vϕ[SUnitRange(joint.minimal_index[2][1], joint.minimal_index[2][2])]
    set_minimal_velocities!(body1, body2, joint, mechanism.timestep, Δv=Δv, Δϕ=Δϕ)
    return body2.state.v15, body2.state.ϕ15
end

################################################################################
# Coordinates and Velocities
################################################################################
function set_minimal_coordinates_velocities!(mechanism::Mechanism, joint::JointConstraint;
        xmin::AbstractVector=szeros(2*control_dimension(joint)))
    pnode = get_body(mechanism, joint.parent_id)
    cnode = get_body(mechanism, joint.child_id)
    set_minimal_coordinates_velocities!(pnode, cnode, joint, mechanism.timestep; xmin=xmin)
end

function set_minimal_coordinates_velocities!(pnode::Node, cnode::Node, joint::JointConstraint, timestep;
        xmin::AbstractVector=szeros(2*control_dimension(joint)))
    nu = control_dimension(joint)
    Δx = xmin[SUnitRange(joint.minimal_index[1]...)]
    Δθ = xmin[SUnitRange(joint.minimal_index[2]...)]
    Δv = xmin[nu .+ SUnitRange(joint.minimal_index[1]...)]
    Δϕ = xmin[nu .+ SUnitRange(joint.minimal_index[2]...)]
	# We need to set the minimal coordinates of the rotational joint first
	# since xb = fct(qb, Δx)
	# since vb = fct(ϕb, Δv)
	set_minimal_coordinates!(pnode, cnode, joint, timestep; Δx=Δx, Δθ=Δθ)
	set_minimal_velocities!(pnode, cnode, joint, timestep; Δv=Δv, Δϕ=Δϕ)
end

function set_minimal_coordinates_velocities_jacobian_minimal(pnode::Node, cnode::Node,
		joint::JointConstraint, timestep)

	xmin = minimal_coordinates_velocities(joint, pnode, cnode, timestep)
	qb = current_orientation(cnode.state)
	function child_maximal_state(xmin)
		set_minimal_coordinates_velocities!(pnode, cnode, joint, timestep; xmin=xmin)
		xb, vb, qb, ϕb = initial_configuration_velocity(cnode.state)
		zb = [xb; vb; vector(qb); ϕb]
		return zb
	end
	J = FiniteDiff.finite_difference_jacobian(xmin -> child_maximal_state(xmin), xmin)
	J = cat(Diagonal(sones(6)), LVᵀmat(qb)', Diagonal(sones(3)), dims=(1,2)) * J
	return J
end

function set_minimal_coordinates_velocities_jacobian_parent(pnode::Node, cnode::Node,
		joint::JointConstraint, timestep)
	xmin = minimal_coordinates_velocities(joint, pnode, cnode, timestep)
	qb = current_orientation(cnode.state)
	xa, va, qa, ϕa = initial_configuration_velocity(pnode.state)
	za = [xa; va; vector(qa); ϕa]

	function child_maximal_state(za)
		xa, va, qa, ϕa = unpack_maximal_state(za,1)
		pnode.state.x2[1] = xa
		pnode.state.v15 = va
		pnode.state.q2[1] = qa
		pnode.state.ϕ15 = ϕa
		set_minimal_coordinates_velocities!(pnode, cnode, joint, timestep; xmin=xmin)
		xb, vb, qb, ϕb = initial_configuration_velocity(cnode.state)
		zb = [xb; vb; vector(qb); ϕb]
		return zb
	end
	J = FiniteDiff.finite_difference_jacobian(za -> child_maximal_state(za), za)
	J = cat(Diagonal(sones(6)), LVᵀmat(qb)', Diagonal(sones(3)), dims=(1,2)) * J
	J = J * cat(Diagonal(sones(6)), LVᵀmat(qa), Diagonal(sones(3)), dims=(1,2))
	return J
end
