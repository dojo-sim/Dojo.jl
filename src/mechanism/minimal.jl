function minimal_to_maximal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, a::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		nu = control_dimension(joint)
		set_minimal_coordinates_velocities!(mechanism,
			joint, xmin=a[off .+ SUnitRange(1,2nu)])
		off += 2nu
	end
	z = get_maximal_state(mechanism)
	return z
end

function minimal_to_maximal_jacobian_analytical(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, a::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	J = zeros(maximal_dimension(mechanism), # TODO: - Nb)
			  minimal_dimension(mechanism))

	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		n = control_dimension(joint)
		idx = collect(off .+ (1:(2n)))
		child_joints = unique([get_node(mechanism, id) for id in recursivedirectchildren!(mechanism.system, joint.id) if get_node(mechanism, id) isa JointConstraint])
		# child_joints = get_child_joints(mechanism, joint)
		function position_velocity(y)
			mech = mechanism

			# currentvals = minimal_coordinates(mech)
			# currentvels = minimal_velocities(mech)
			cv = minimal_coordinates_velocities(mech)

			# set_joint_position!(mech, joint, y[1:n])
			# set_velocity!(mech, joint, y[n .+ (1:n)])
			set_minimal_coordinates_velocities!(mechanism, joint, xmin=y[off .+ SUnitRange(1,2n)])

			for node in child_joints
				# set_joint_position!(mech, node, currentvals[node.id])
				# set_velocity!(mech, node, currentvels[node.id])
				set_minimal_coordinates_velocities!(mechanism, joint, xmin=cv[node.id])
			end

			return get_maximal_state(mech)
		end

		J[:, idx] = FiniteDiff.finite_difference_jacobian(position_velocity, a[idx])

		function joint_position_velocity(mech, joint, θ)
			n = control_dimension(joint)
			# x, q = set_joint_position!(mech, joint, θ[1:n])
			# v, ω = set_velocity!(mech, joint, θ[n .+ (1:n)])
			set_minimal_coordinates_velocities!(mechanism, joint, xmin=θ[SUnitRange(1,2n)])
			x2, v15, q2, ϕ15 = initial_configuration_velocity(get_body(joint.child_id).state)
			return [x2; v15; vector(q2); ϕ15]
		end

		function joint_position_velocity(mech, joint, z, θ)
			n = control_dimension(joint)

			body_parent = get_body(mech, joint.parent_id)
			xp = z[1:3]
			vp = z[4:6]
			qp = UnitQuaternion(z[7:10]..., false)
			ϕp = z[11:13]

			set_position!(body_parent, x=xp, q=qp)
			set_velocity!(body_parent, v=vp, ω=ϕp)

			x, q = set_joint_position!(mech, joint, θ[1:n])
			v, ω = set_velocity!(mech, joint, θ[n .+ (1:n)])

			return [x; v; vector(q); ω]
		end

		function joint_position_velocity_jacobian(mech, joint, θ)
			FiniteDiff.finite_difference_jacobian(y -> joint_position_velocity(mech, joint, y), θ)
		end

		function position_velocity_jacobian(θ)
			G = zeros(maximal_dimension(mechanism), length(θ))

			mech = mechanism

			currentvals = minimal_coordinates(mech)
			currentvels = minimal_velocities(mech)

			x, q = set_joint_position!(mech, joint, θ[1:n])
			v, ω = set_velocity!(mech, joint, θ[n .+ (1:n)])
			zp = [x; v; vector(q); ω]

			for node in child_joints
				set_joint_position!(mech, node, currentvals[node.id])
				set_velocity!(mech, node, currentvels[node.id])
			end

			# root
			∂z∂θ = joint_position_velocity_jacobian(mech, joint, θ)
			G[(joint.child_id - 1 - Ne) * 13 .+ (1:13), :] = ∂z∂θ

			# branch
			y = ∂z∂θ
			for node in child_joints
				∂z∂z = FiniteDiff.finite_difference_jacobian(b -> joint_position_velocity(mech, node, b, [currentvals[node.id]; currentvels[node.id]]), zp)
				y = ∂z∂z * ∂z∂θ
				G[(node.child_id - 1 - Ne) * 13 .+ (1:13), :] = y
				zp = joint_position_velocity(mech, node, zp, θ)
			end

			return G
		end

		# J[:, idx] = position_velocity_jacobian(a[idx])

		off += 2n
	end

	return J
end

function minimal_to_maximal_jacobian(mechanism::Mechanism, x)
	FiniteDiff.finite_difference_jacobian(y -> minimal_to_maximal(mechanism, y), x)
end

function get_minimal_gradients(mechanism::Mechanism{T}, z::AbstractVector{T}, u::AbstractVector{T};
	opts=SolverOptions()) where T
	# simulate next state
	step!(mechanism, z, u, opts=opts)
	# current maximal state
	z = get_state(mechanism)
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

function get_minimal_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni};
	pos_noise=nothing, vel_noise=nothing,
	pos_noise_range=[-Inf, Inf], vel_noise_range=[-3.9 / mechanism.timestep^2, 3.9 / mechanism.timestep^2]) where {T,Nn,Ne,Nb,Ni}
	x = []
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		c = zeros(T,0)
		v = zeros(T,0)
		pbody = get_body(mechanism, joint.parent_id)
		cbody = get_body(mechanism, joint.child_id)
		for (i, element) in enumerate(joint.constraints)
			pos = minimal_coordinates(element, pbody, cbody)
			vel = minimal_velocities(element, pbody, cbody)
			if pos_noise != nothing
				pos += clamp.(length(pos) == 1 ? rand(pos_noise, length(pos))[1] : rand(pos_noise, length(pos)), pos_noise_range...)
			end
			if vel_noise != nothing
				vel += clamp.(length(vel) == 1 ? rand(vel_noise, length(vel))[1] : rand(vel_noise, length(vel)), vel_noise_range...)
			end
			push!(c, pos...)
			push!(v, vel...)
		end
		push!(x, [c; v]...)
	end
	x = [x...]
	return x
end
