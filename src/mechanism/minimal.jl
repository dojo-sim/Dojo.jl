function minimal_to_maximal(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector{Tx}) where {T,Nn,Ne,Nb,Ni,Tx}
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	off = 0
	for id in reverse(mechanism.system.dfs_list)
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
		n = control_dimension(joint)
		if joint.parent_id != nothing
			c = x[off .+ (1:n)]; off += n
			v = x[off .+ (1:n)]; off += n # in body1
			set_joint_position!(mechanism, joint, c)
			set_velocity!(mechanism, joint, v) # in body1
		else
			@assert length(Set(joint.child_ids)) == 1 # only one body is linked to the origin
			c = zeros(Tx,0)
			v = zeros(Tx,0)
			# we need a special case: when the first link has free rotation wrt the origin
			q2 = one(UnitQuaternion) # only use for the special case
			for element in joint.constraints
				nj = control_dimension(element)
				push!(c, x[off .+ (1:nj)]...); off += nj
			end
			push!(v, x[off .+ (1:n)]...); off += n
			set_joint_position!(mechanism, joint, c)
			set_velocity!(mechanism, joint, v) # assumes we provide v and ϕ in body1's coordinates i.e world coordinates
		end
	end
	z = get_maximal_state(mechanism)
	return z
end

function minimal_to_maximal_jacobian(mechanism::Mechanism, x)
	FiniteDiff.finite_difference_jacobian(x -> minimal_to_maximal(mechanism, x), x)
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
		for (i,element) in enumerate(joint.constraints)
			cbody = get_body(mechanism, joint.child_ids[i])
			pbody = get_body(mechanism, joint.parent_id)
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