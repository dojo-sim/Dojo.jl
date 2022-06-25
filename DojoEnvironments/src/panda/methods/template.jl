	# Inverse kinematics: will try to find a configuration that set the first contact point at p_end_effector
	# at the p_foot location
	function panda_inverse_kinematics(mechanism::Mechanism, p_end_effector, q_end_effector;
			θ0=[π, -0.8, 0.0, 1.6, 0.0, -2.4, 0.0])
		# starting point of the local search

		nu = input_dimension(mechanism)
		θ = [θ0; zeros(nu - 7)]

		for k = 1:100
			err = panda_inverse_kinematics_error(mechanism, p_end_effector, q_end_effector, θ)
			norm(err, Inf) < 1e-10 && continue
			∇ = FiniteDiff.finite_difference_jacobian(θ -> panda_inverse_kinematics_error(mechanism, p_end_effector, q_end_effector, θ), θ)
			θ -= (∇' * inv(∇*∇' + 1e-2*I)) * err
		end
		x = vcat([[θi, 0.0] for θi in θ]...)
		z = minimal_to_maximal(mechanism, x)
		return z, θ[1:7]
	end

	function panda_inverse_kinematics_error(mechanism::Mechanism, p_end_effector, q_end_effector, θ)

		set_minimal_state!(mechanism, vcat([[θi, 0.0] for θi in θ]...))
		contact = mechanism.contacts[1]
		body = get_body(mechanism, contact.parent_id)
		p = contact_location(mechanism, contact)
		q = current_orientation(body.state)

		err = [p - p_end_effector; vector(q) - vector(q_end_effector)]

		return err
	end

	function panda_inverse_kinematics_trajectory(mechanism::Mechanism{T},
			p_end_effector::Vector{SVector{3,T}}, q_end_effector::Vector{Quaternion{T}};
			θ::Vector{T}=[π, -0.8, 0.0, 1.6, 0.0, -2.4, 0.0]) where T
		traj = []
		N = length(p_end_effector)
		for i = 1:N
			@show i, N
			z, θ = panda_inverse_kinematics(mechanism::Mechanism, p_end_effector[i], q_end_effector[i]; θ0=θ)
			push!(traj, z)
		end
		return traj
	end
