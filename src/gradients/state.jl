"""
	maximal_to_minimal_jacobian(mechanism, z)

	Jacobian of mapping from maximal to minimal representation

	mechanism: Mechanism
	z: maximal state
"""
function maximal_to_minimal_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	J = zeros(minimal_dimension(mechanism), maximal_dimension(mechanism) - Nb)
	timestep= mechanism.timestep
	row_shift = 0
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c_shift = 0
		v_shift = input_dimension(joint)
		ichild = joint.child_id - Ne
		for element in [joint.translational, joint.rotational]
			nu_element = input_dimension(element)

			c_idx = row_shift + c_shift .+ (1:nu_element)
			v_idx = row_shift + v_shift .+ (1:nu_element)

			xb, vb, qb, ϕb = unpack_maximal_state(z, ichild)

			xb_idx = collect((ichild-1)*12 .+ (1:3))
			vb_idx = collect((ichild-1)*12 .+ (4:6))
			qb_idx = collect((ichild-1)*12 .+ (7:9))
			ϕb_idx = collect((ichild-1)*12 .+ (10:12))

			if joint.parent_id != 0
				iparent = joint.parent_id - Ne
				xa, va, qa, ϕa = unpack_maximal_state(z, iparent)

				xa_idx = collect((iparent-1)*12 .+ (1:3))
				va_idx = collect((iparent-1)*12 .+ (4:6))
				qa_idx = collect((iparent-1)*12 .+ (7:9))
				ϕa_idx = collect((iparent-1)*12 .+ (10:12))

				J[c_idx, [xa_idx; qa_idx]] = minimal_coordinates_jacobian_configuration(:parent, element, xa, qa, xb, qb)
				J[v_idx, [xa_idx; qa_idx]] = minimal_velocities_jacobian_configuration(:parent, element, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
				J[v_idx, [va_idx; ϕa_idx]] = minimal_velocities_jacobian_velocity(:parent, element, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
			else
				xa, va, qa, ϕa = current_configuration_velocity(mechanism.origin.state)
			end

			J[c_idx, [xb_idx; qb_idx]] = minimal_coordinates_jacobian_configuration(:child, element, xa, qa, xb, qb)
			J[v_idx, [xb_idx; qb_idx]] = minimal_velocities_jacobian_configuration(:child, element, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)
			J[v_idx, [vb_idx; ϕb_idx]] = minimal_velocities_jacobian_velocity(:child, element, xa, va, qa, ϕa, xb, vb, qb, ϕb, timestep)

			c_shift += nu_element
			v_shift += nu_element
		end
		row_shift += 2 * input_dimension(joint)
	end
	return J
end

"""
    get_maximal_gradients!(mechanism, z, u; opts)

    return maximal gradients for mechanism
    note: this requires simulating the mechanism for one time step

    mechanism: Mechanism
    z: state
    u: input
    opts: SolverOptions
"""
function get_maximal_gradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
    opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}

    step!(mechanism, z, u, opts=opts)
    jacobian_state, jacobian_control = get_maximal_gradients(mechanism)

    return jacobian_state, jacobian_control
end

function get_maximal_gradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	nu = input_dimension(mechanism)

	for entry in mechanism.data_matrix.nzval # reset matrix
		entry.value .= 0.0
	end
	jacobian_data!(mechanism.data_matrix, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = data_dim.(nodes)
	index_row = [1+sum(dimrow[1:i-1]):sum(dimrow[1:i]) for i in 1:length(dimrow)]
	index_col = [1+sum(dimcol[1:i-1]):sum(dimcol[1:i]) for i in 1:length(dimcol)]

	index_state = [index_col[body.id][[14:16; 8:10; 17:19; 11:13]] for body in mechanism.bodies] # ∂ x2 v15 q2 ϕ15
	index_control = [index_col[joint.id][1:input_dimension(joint)] for joint in mechanism.joints] # ∂ u

	datamat = full_matrix(mechanism.data_matrix, dimrow, dimcol)
	solmat = full_matrix(mechanism.system)

	# data Jacobian
	data_jacobian = solmat \ datamat #TODO: use pre-factorization

	# Jacobian
	jacobian_state = zeros(12Nb,12Nb)
	jacobian_control = zeros(12Nb,nu)
	for (i, body) in enumerate(mechanism.bodies)
		id = body.id
		# Fill in gradients of v25, ϕ25
		jacobian_state[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_state...)]
		jacobian_control[12*(i-1) .+ [4:6; 10:12],:] += data_jacobian[index_row[id], vcat(index_control...)]

		# Fill in gradients of x3, q3
		x2 = body.state.x2
		q2 = body.state.q2
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		q3 = next_orientation(q2, ϕ25, timestep)
		jacobian_state[12*(i-1) .+ (1:3), :] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (1:3), 12*(i-1) .+ (1:3)] += linear_integrator_jacobian_position(x2, v25, timestep)
		jacobian_state[12*(i-1) .+ (7:9), :] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_state...)]
		jacobian_state[12*(i-1) .+ (7:9), 12*(i-1) .+ (7:9)] += LVᵀmat(q3)' * rotational_integrator_jacobian_orientation(q2, ϕ25, timestep, attjac=true)

		jacobian_control[12*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(x2, v25, timestep) * data_jacobian[index_row[id][1:3], vcat(index_control...)]
		jacobian_control[12*(i-1) .+ (7:9),:] += LVᵀmat(q3)' * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian[index_row[id][4:6], vcat(index_control...)]
	end

	return jacobian_state, jacobian_control
end

"""
	minimal_to_maximal_jacobian(mechanism, x)

	Jacobian of mapping from minimal to maximal representation

	mechanism: Mechanism
	y: minimal state
"""
function minimal_to_maximal_jacobian(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	timestep= mechanism.timestep
	J = zeros(maximal_dimension(mechanism, attjac=true), minimal_dimension(mechanism))

	# Compute partials
	partials = Dict{Vector{Int}, Matrix{T}}()
	for cnode in mechanism.bodies
		for joint in parent_joints(mechanism, cnode)
			pnode = get_node(mechanism, joint.parent_id, origin=true)
			partials[[cnode.id, joint.id]] = set_minimal_coordinates_velocities_jacobian_minimal(joint, pnode, cnode, timestep) # 12 x 2nu (xvqϕ x Δxθvϕ)
			partials[[cnode.id, pnode.id]] = set_minimal_coordinates_velocities_jacobian_parent(joint, pnode, cnode, timestep) # 12 x 12 (xvqϕ x xvqϕ)
		end
	end

	# Index
	row = [12(i-1)+1:12i for i = 1:Nb]
	col = [] # ordering joints from root to tree
	col_idx = zeros(Int,Ne)
	cnt = 0
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only keep joints
		cnt += 1
		nu = input_dimension(get_joint(mechanism, id))
		if length(col) > 0
			push!(col, col[end][end] .+ (1:2nu))
		else
			push!(col, 1:2nu)
		end
		col_idx[id] = cnt
	end

	 # chain partials together from root to leaves
	for id in mechanism.root_to_leaves
		!(Ne < id <= Ne+Nb) && continue # only treat bodies
		cnode = get_node(mechanism, id)
		for joint in parent_joints(mechanism, cnode)
			pnode = get_node(mechanism, joint.parent_id, origin=true)
			J[row[cnode.id-Ne], col[col_idx[joint.id]]] += partials[[cnode.id, joint.id]] # ∂zi∂θp(i)
			(pnode.id == 0) && continue # avoid origin
			J[row[cnode.id-Ne], :] += partials[[cnode.id, pnode.id]] * J[row[pnode.id-Ne], :] # ∂zi∂zp(p(i)) * ∂zp(p(i))/∂θ
		end
	end
	return J
end

"""
    get_minimal_gradients!(mechanism, y, u; opts)

    return minimal gradients for mechanism
    note: this requires simulating the mechanism for one time step

    mechanism: Mechanism
    y: state
    u: input
	opts: SolverOptions
"""
function get_minimal_gradients!(mechanism::Mechanism{T}, y::AbstractVector{T}, u::AbstractVector{T};
	opts=SolverOptions()) where T
	# simulate next state
	step!(mechanism, y, u, opts=opts)
	return get_minimal_gradients!(mechanism)
end

function get_minimal_gradients!(mechanism::Mechanism{T}) where T
	# current maximal state
	z = get_maximal_state(mechanism)
	# next maximal state
	z_next = get_next_state(mechanism)
	# current minimal state
	x = maximal_to_minimal(mechanism, z)
	# maximal dynamics Jacobians
	maximal_jacobian_state, minimal_jacobian_control = get_maximal_gradients(mechanism)
	# minimal to maximal Jacobian at current time step (rhs)
	min_to_max_jacobian_current = minimal_to_maximal_jacobian(mechanism, x)
	# maximal to minimal Jacobian at next time step (lhs)
	max_to_min_jacobian_next = maximal_to_minimal_jacobian(mechanism, z_next)
	# minimal state Jacobian
	minimal_jacobian_state = max_to_min_jacobian_next * maximal_jacobian_state * min_to_max_jacobian_current
	# minimal control Jacobian
	minimal_jacobian_control = max_to_min_jacobian_next * minimal_jacobian_control
	return minimal_jacobian_state, minimal_jacobian_control
end
