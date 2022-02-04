function get_maximal_gradients(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep = mechanism.timestep
	nu = control_dimension(mechanism)
	attjac = false
	nic = attjac ? 12Nb : 13Nb
	njoints = joint_dimension(mechanism)
	datamat = full_data_matrix(mechanism, attjac = attjac)
	solmat = full_matrix(mechanism.system)

	# data Jacobian 
	data_jacobian = - solmat \ datamat #TODO: use pre-factorization
	data_jacobian_state = data_jacobian[njoints .+ (1:6Nb),1:nic]
	data_jacobian_control = data_jacobian[njoints .+ (1:6Nb),nic .+ (1:nu)]

	# Jacobian
	jacobian_state = zeros(13Nb,13Nb)
	jacobian_control = zeros(13Nb,nu)
	for (i, body) in enumerate(mechanism.bodies)
		# Fill in gradients of v25, ϕ25
		jacobian_state[13*(i-1) .+ [4:6; 11:13],:] += data_jacobian_state[6*(i-1) .+ (1:6),:]
		jacobian_control[13*(i-1) .+ [4:6; 11:13],:] += data_jacobian_control[6*(i-1) .+ (1:6),:]

		# Fill in gradients of x3, q3
		q2 = body.state.q2[1]
		ϕ25 = body.state.ϕsol[2]
		jacobian_state[13*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(timestep) * data_jacobian_state[6*(i-1) .+ (1:3),:]
		jacobian_state[13*(i-1) .+ (1:3),13*(i-1) .+ (1:13)] += linear_integrator_jacobian_position() * [I(3) zeros(3,10)]
		jacobian_state[13*(i-1) .+ (7:10),:] += rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian_state[6*(i-1) .+ (4:6),:]
		jacobian_state[13*(i-1) .+ (7:10),13*(i-1) .+ (1:13)] += rotational_integrator_jacobian_orientation(q2, ϕ25, timestep, attjac = false) * [zeros(4,6) I(4) zeros(4,3)]

		jacobian_control[13*(i-1) .+ (1:3),:] += linear_integrator_jacobian_velocity(timestep) * data_jacobian_control[6*(i-1) .+ (1:3),:]
		jacobian_control[13*(i-1) .+ (7:10),:] += rotational_integrator_jacobian_velocity(q2, ϕ25, timestep) * data_jacobian_control[6*(i-1) .+ (4:6),:]
	end
	return jacobian_state, jacobian_control
end

function get_maximal_gradients!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{T}, u::AbstractVector{T};
    opts=SolverOptions()) where {T,Nn,Ne,Nb,Ni}
    step!(mechanism, z, u, opts=opts)
    jacobian_state, jacobian_control = get_maximal_gradients(mechanism)
    return jacobian_state, jacobian_control
end

function maximal_to_minimal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	x = []
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c = zeros(Tz,0)
		v = zeros(Tz,0)
		ichild = joint.child_id - Ne
		for (i,element) in enumerate(joint.constraints)
			xb, vb, qb, ϕb = unpack_maximal_state(z, ichild)
			if joint.parent_id != 0
				iparent = joint.parent_id - Ne
				xa, va, qa, ϕa = unpack_maximal_state(z, iparent)
			else
				xa, va, qa, ϕa = current_configuration_velocity(mechanism.origin.state)
			end
			if typeof(element) <: Translational
				push!(c, minimal_coordinates(element, xa, qa, xb, qb)...) # Δx in bodya's coordinates projected on elementAB's nullspace
				push!(v, minimal_velocities(element, xa, va, qa, ϕa, xb, vb, qb, ϕb)...) # Δv in bodya's coordinates projected on elementAB's nullspace
			elseif typeof(element) <: Rotational
				push!(c, minimal_coordinates(element, qa, qb)...) # Δq in bodya's coordinates projected on elementAB's nullspace
				push!(v, minimal_velocities(element, qa, ϕa, qb, ϕb)...) # Δϕ in bodya's coordinates projected on jointAB's nullspace
			end
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end

function maximal_to_minimal_jacobian_analytical(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, z::AbstractVector{Tz}) where {T,Nn,Ne,Nb,Ni,Tz}
	J = zeros(minimal_dimension(mechanism), maximal_dimension(mechanism))# TODO: - Nb)
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	row_shift = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c_shift = 0 
		v_shift = control_dimension(joint)
		ichild = joint.child_id - Ne
		for (i, element) in enumerate(joint.constraints)

			xb, vb, qb, ϕb = unpack_maximal_state(z, ichild)

			# xb_idx = collect((ichild-1)*12 .+ (1:3)) 
			# vb_idx = collect((ichild-1)*12 .+ (4:6)) 
			# qb_idx = collect((ichild-1)*12 .+ (7:9)) 
			# ϕb_idx = collect((ichild-1)*12 .+ (10:12))
			xb_idx = collect((ichild-1)*13 .+ (1:3)) 
			vb_idx = collect((ichild-1)*13 .+ (4:6)) 
			qb_idx = collect((ichild-1)*13 .+ (7:10)) 
			ϕb_idx = collect((ichild-1)*13 .+ (11:13))

			if joint.parent_id != 0
				iparent = joint.parent_id - Ne
				xa, va, qa, ϕa = unpack_maximal_state(z, iparent)

				# xa_idx = collect((iparent-1)*12 .+ (1:3)) 
				# va_idx = collect((iparent-1)*12 .+ (4:6)) 
				# qa_idx = collect((iparent-1)*12 .+ (7:9)) 
				# ϕa_idx = collect((iparent-1)*12 .+ (10:12)) 
				xa_idx = collect((iparent-1)*13 .+ (1:3)) 
				va_idx = collect((iparent-1)*13 .+ (4:6)) 
				qa_idx = collect((iparent-1)*13 .+ (7:10)) 
				ϕa_idx = collect((iparent-1)*13 .+ (11:13)) 
			else
				xa, va, qa, ϕa = current_configuration_velocity(mechanism.origin.state)
			end 

			nu_element = control_dimension(element)

			c_idx = row_shift + c_shift .+ (1:nu_element)
			v_idx = row_shift + v_shift .+ (1:nu_element)

			if typeof(element) <: Translational 
				if joint.parent_id != 0
					J[c_idx, [xa_idx; qa_idx]] = minimal_coordinates_jacobian_configuration(:parent, element, xa, qa, xb, qb)
					J[v_idx, [xa_idx; qa_idx]] = minimal_velocities_jacobian_configuration(:parent, element, xa, va, qa, ϕa, xb, vb, qb, ϕb)
					J[v_idx, [va_idx; ϕa_idx]] = minimal_velocities_jacobian_velocity(:parent, element, xa, va, qa, ϕa, xb, vb, qb, ϕb)
				end

				J[c_idx, [xb_idx; qb_idx]] = minimal_coordinates_jacobian_configuration(:child, element, xa, qa, xb, qb)
				J[v_idx, [xb_idx; qb_idx]] = minimal_velocities_jacobian_configuration(:child, element, xa, va, qa, ϕa, xb, vb, qb, ϕb)
				J[v_idx, [vb_idx; ϕb_idx]] = minimal_velocities_jacobian_velocity(:child, element, xa, va, qa, ϕa, xb, vb, qb, ϕb)

			elseif typeof(element) <: Rotational
				if joint.parent_id != 0
					J[c_idx, [xa_idx; qa_idx]] = minimal_coordinates_jacobian_configuration(:parent, element, qa, qb)
					J[v_idx, [xa_idx; qa_idx]] = minimal_velocities_jacobian_configuration(:parent, element, qa, ϕa, qb, ϕb)
					J[v_idx, [va_idx; ϕa_idx]] = minimal_velocities_jacobian_velocity(:parent, element, qa, ϕa, qb, ϕb)
				end
				J[c_idx, [xb_idx; qb_idx]] = minimal_coordinates_jacobian_configuration(:child, element, qa, qb)
				J[v_idx, [xb_idx; qb_idx]] = minimal_velocities_jacobian_configuration(:child, element, qa, ϕa, qb, ϕb)
				J[v_idx, [vb_idx; ϕb_idx]] = minimal_velocities_jacobian_velocity(:child, element, qa, ϕa, qb, ϕb)

			end
			c_shift += nu_element
			v_shift += nu_element
		end
		row_shift += 2 * control_dimension(joint)
	end
	return J
end

function maximal_to_minimal_jacobian(mechanism::Mechanism, z)
	FiniteDiff.finite_difference_jacobian(y -> maximal_to_minimal(mechanism, y), z)
end

function get_maximal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T, 13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		x2 = body.state.x2[1]
		v15 = body.state.v15
		q2 = body.state.q2[1]
		ϕ15 = body.state.ϕ15
		set_maximal_state!(z, x2, v15, q2, ϕ15, i)
	end
	return z
end

function unpack_maximal_state(z::AbstractVector, i::Int)
	zi = z[(i-1)*13 .+ (1:13)]
	x2 = zi[1:3]
	v15 = zi[4:6]
	q2 = UnitQuaternion(zi[7:10]..., false)
	ϕ15 = zi[11:13]
	return x2, v15, q2, ϕ15
end

function set_maximal_state!(z::AbstractVector, x2::AbstractVector, v15::AbstractVector,
		q2::UnitQuaternion, ϕ15::AbstractVector, i::Int)
	z[(i-1)*13 .+ (1:3)] = x2
	z[(i-1)*13 .+ (4:6)] = v15
	z[(i-1)*13 .+ (7:10)] = vector(q2)
	z[(i-1)*13 .+ (11:13)] = ϕ15
	return nothing
end

