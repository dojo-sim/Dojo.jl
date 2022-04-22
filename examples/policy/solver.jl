function augmented_lagrangian_solver(env, obj, Xref, splits, vars)
	N = length(splits)
	x0 = Xref[1]
	X = policy_rollout_only(env, x0, splits, vars)
	for i = 1:20
		(violation(X, N) < 1e-4) && break
		println("")
		println("i j  f       g       α       Δ       vio")
		for j = 1:4
			f = cost_only(env, obj, Xref, splits, vars)
			g = FiniteDiff.finite_difference_gradient(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			# g = ForwardDiff.gradient(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			(norm(g, Inf) < 1e-3) && break
			hess = FiniteDiff.finite_difference_hessian(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			# hess = ForwardDiff.hessian(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			Δ = - hess \ g
			X = policy_rollout_only(env, x0, splits, vars)
			vio = violation(X, N)
			α = linesearch(env, obj, Xref, splits, vars, Δ, f, vio)
			vars += α * Δ
			println("$i $j " *
				"$(scn(f, exp_digits=2))" *
				"$(scn(norm(g, Inf), exp_digits=2))" *
				"$(scn(α, exp_digits=2))" *
				"$(scn(norm(Δ, Inf), exp_digits=2))" *
				"$(scn(vio, exp_digits=2))"
				)
		end
		# Dual ascent and penalty update
		X = policy_rollout_only(env, x0, splits, vars)
		plot_rollout(X, Xref, splits)
		con = constraints(X, N)
		for i = 1:N-1
			obj.λ[i] += obj.ρ * con[i]
		end
		obj.ρ = min(obj.ρ * 10, 1e6)
	end
	return vars
end

function linesearch(env, obj, Xref, splits, vars, Δ, f, vio)
	N = length(splits)
	α = 1.0
	for i = 1:10
		f1 = cost_only(env, obj, Xref, splits, vars+α*Δ)
		X1 = policy_rollout_only(env, x0, splits, vars+α*Δ)
		(f1 <= f || violation(X1, N) <= vio) && break
		α *= 0.5
	end
	return α
end

function continuation_solve(env, obj, Xref, vars, H, N)
	x0 = Xref[1]
	obj.ρ = 1e-3
	for Ni in [20, 10, 5, 2]
		splits = split_trajectory(H, Ni)
		obj.λ = [zeros(nx) for i=1:Ni-1]
		obj.Q *= 3

		θ, xstarts = unpack_vars(vars, N=Ni)
		vars = [θ; vcat([Xref[s[end]] for s in splits[1:Ni-1]]...)]
		vars = augmented_lagrangian_solver(env, obj, Xref, splits, vars)
	end
	return vars
end
