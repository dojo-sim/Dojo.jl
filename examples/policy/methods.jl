function policy_rollout_only(env::Environment, x0, splits, vars::Vector)
	N = length(splits)
	H = splits[end][end]
	θ, xstarts,us = unpack_vars(vars, N=N, H=H)

	X = [[zeros(nx) for j=1:length(splits[i])] for i=1:N]

	# Initial states
	X[1][1] .= deepcopy(x0)
	for i = 1:N-1
		X[i+1][1] .= deepcopy(xstarts[i])
	end
	for i = 1:N
		for (j,ind) in enumerate(splits[i][1:end-1])
			x = X[i][j]
			y = X[i][j+1]
			xref = Xref[ind]
			# u = policy(env, x, θ)
			u = us[ind]
			# dynamics(y, env, x, u, nothing;
			# 	attitude_decompress=false)
			dx = zeros(nx,nx)
			du = zeros(nx,nu)
			dummy_dynamics(y, dx, du, x, u)
		end
	end
    return X
end

function pure_policy_rollout(env::Environment, x0, splits, vars::Vector; mode::Symbol=:policy)
	N = length(splits)
	H = splits[end][end]
	θ, xstarts, us = unpack_vars(vars, N=N, H=H)

	X = [zeros(nx) for i = 1:H]
	X[1] .= deepcopy(x0)
	for i = 1:H-1
		x = X[i]
		y = X[i+1]
		if mode == :policy
			u = policy(env, x, θ)
		else mode == :control
			u = us[i]
		end
		dx = zeros(nx,nx)
		du = zeros(nx,nu)
		dummy_dynamics(y, dx, du, x, u)
	end
    return X
end

function cost_only(env::Environment, obj, Xref, splits, vars)
	X = policy_rollout_only(env, Xref[1], splits, vars)
	cost_only(env, obj, X, Xref, splits, vars)
end


function cost_only(env::Environment, obj, X, Xref, splits, vars)
	N = length(splits)
	H = splits[end][end]
	θ, xstarts, us = unpack_vars(vars, N=N, H=H)

	# Initialization
	f = 0.0

	# Evaluation
	# Objective
	for i = 1:N
		for (j,ind) in enumerate(splits[i])
			x = X[i][j]
			xref = Xref[ind]
			f += 0.5 * (x - xref)' * obj.Q * (x - xref)
			if ind < H
				up = policy(env, x, θ)
				uc = us[ind]
				f += obj.λu[ind]' * (up - uc)
				f += 0.5 * obj.ρu * (up - uc)' * (up - uc)
				f += 0.5 * up' * obj.R * up
			end
		end
	end

	# Constraints
	for i = 1:N-1
		λx = obj.λx[i]
		xend = X[i][end]
		xstart = X[i+1][1]
		f += λx' * (xstart - xend)
		f += 0.5 * obj.ρx * (xstart - xend)' * (xstart - xend)
		# xref = Xref[splits[i][end]]
		# f += λx' * (xref - xend)
		# f += 0.5 * obj.ρ * (xref - xend)' * (xref - xend)
	end

	# parameter regularizer
	for i = 1:N-1
		xstart = xstarts[i]
		xref = Xref[splits[i][end]]
		f += 0.5 * (xstart - xref)' * obj.Qstart * (xstart - xref)
		# f += 0.5 * (xstart - Xref[1])' * obj.Qstart * (xstart - Xref[1])
	end
	f += 0.5 * θ' * obj.Qθ * θ

	return f
end

################################################################################
# vars0 = zeros(nva)
# θ, xs = unpack_vars(vars0)
# X0, dX0, dU0 = policy_rollout(env, x0, vars0)
#
# Q0 = 0*1e-0*Diagonal(ones(nx))
# R0 = 0*1e-2*Diagonal(ones(nu))
# λ0 = [(i in (1,2))*ones(nx) for i=1:N-1]
# ρ0 = 0.0
# obj0 = AugmentedObjective119(Q0, R0, λ0, ρ0)
#
# cost(env, obj0, X0, dX0, dU0, Xref, vars0)
# cost(env, obj0, Xref, vars0)
#
# g1 = FiniteDiff.finite_difference_gradient(vars -> cost_only(env, obj0, Xref, vars), vars0)
# f0, g0 = cost(env, obj0, Xref, vars0)
# norm(g0 - g1)







function policy_rollout(env::Environment, x0, vars::Vector)
	θ, xstarts = unpack_vars(vars)
	splits = split_trajectory(H, N)

	X = [[zeros(nx) for j=1:length(splits[i])] for i=1:N]
	dX = [[zeros(nx,nx) for j=1:length(splits[i])] for i=1:N]
	dU = [[zeros(nx,nu) for j=1:length(splits[i])] for i=1:N]

	# Initial states
	X[1][1] .= deepcopy(x0)
	for i = 1:N-1
		X[i+1][1] .= deepcopy(xstarts[i])
	end
	for i = 1:N
		for (j,ind) in enumerate(splits[i][1:end-1])
			x = X[i][j]
			y = X[i][j+1]
			dx = dX[i][j+1]
			du = dU[i][j+1]
			xref = Xref[ind]
			u = policy(env, x, θ)
			dummy_dynamics(y, dx, du, x, u)
			# dynamics_and_jacobian(y, dx, du, env, x, u, nothing;
				# attitude_decompress=false)
		end
	end
    return X, dX, dU
end

function cost(env::Environment, obj, Xref, vars)
	X, dX, dU = policy_rollout(env, Xref[1], vars)
	cost(env, obj, X, dX, dU, Xref, vars)
end

function cost(env::Environment, obj, X, dX, dU, Xref, vars)
	splits = split_trajectory(H, N)
	θ, xstarts = unpack_vars(vars, N=N)

	# Initialization
	f = 0.0
	g_θ = zeros(nθ)
	g_xstarts = [zeros(nx) for i=1:N-1]

	## Evaluation
	# Objective
	for i = 1:N
		for (j,ind) in enumerate(splits[i])
			x = X[i][j]
			xref = Xref[ind]
			f += 0.5 * (x - xref)' * obj.Q * (x - xref)
		end
	end
	# Constraints
	for i = 1:N-1
		λ = obj.λ[i]
		xend = X[i][end]
		xstart = X[i+1][1]
		f += λ' * (xstart - xend)
		f += 0.5 * obj.ρ * (xstart - xend)' * (xstart - xend)
	end

	## Gradient
	for i = 1:N-1
		λ = obj.λ[i]
		ind = splits[i+1][1] # starting index
		xref = Xref[ind]
		xend = X[i][end]
		xstart = X[i+1][1]
		# Objective
		g_xstarts[i] += obj.Q * (xstart - xref)
		# Constraints
		g_xstarts[i] += λ + obj.ρ * (xstart - xend)
	end
	for i = 1:N
		∂x∂θ = zeros(nx, nθ)
		for (j,ind) in enumerate(splits[i][1:end-1])
			x0 = X[i][j] # x0
			x1 = X[i][j+1] # x1
			x1ref = Xref[ind+1]

			u = policy(env, x0, θ)
			∂u∂θ = policy_jacobian_parameters(env, x0, θ)
			∂u∂x0 = policy_jacobian_state(env, x0, θ)

			∂x1∂u = dU[i][j+1]
			∂x1∂x0 = dX[i][j+1]

			∂x1∂θ = ∂x1∂u * ∂u∂θ
			∂x1∂θ += (∂x1∂x0 + ∂x1∂u * ∂u∂x0) * ∂x∂θ

			∂o∂x1 = (x1 - x1ref)' * obj.Q
			∂o∂u = u' * obj.R * u
			g_θ += (∂o∂x1 * ∂x1∂θ + ∂o∂u * ∂u∂θ)'

			# For the last element we have a stitching constraint, except for the last split
			if j == length(splits[i])-1 && i < N
				xend = x1
				xstart = X[i+1][1]
				λ = obj.λ[i]
				∂xend∂θ = ∂x1∂θ
				∂c∂xend = -λ' + obj.ρ * (xend - xstart)' # xend = x1
				g_θ += (∂c∂xend * ∂xend∂θ)'
			end
			# update
			∂x∂θ = ∂x1∂θ
		end
	end
	g_θ = g_θ[:,1]
	g = pack_vars(g_θ, g_xstarts, N=N)
	return f, g
end
