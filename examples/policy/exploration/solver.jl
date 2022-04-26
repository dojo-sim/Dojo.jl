function augmented_lagrangian_solver(env, obj, Xref, splits, vars)
	N = length(splits)
	x0 = Xref[1]
	X = policy_rollout_only(env, x0, splits, vars)
	for i = 1:20
		(violation(X, obj, splits, vars) < 1e-3) && break
		println("")
		println("i j  f       g       α       Δ       vio")
		for j = 1:10
			f = cost_only(env, obj, Xref, splits, vars)
			g = FiniteDiff.finite_difference_gradient(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			# g = ForwardDiff.gradient(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			(norm(g, Inf) < 1e-4) && break
			hess = FiniteDiff.finite_difference_hessian(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			# hess = ForwardDiff.hessian(vars -> cost_only(env, obj, Xref, splits, vars), vars)
			Δ = - hess \ g
			X = policy_rollout_only(env, x0, splits, vars)
			vio = violation(X, obj, splits, vars)
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
		dual_ascent(env, obj, Xref, x0, splits, vars)
		obj.ρx = min(obj.ρx * 3, 1e6)
		obj.ρu = min(obj.ρu * 3, 1e6)
	end
	return vars
end

function dual_ascent(env, obj, Xref, x0, splits, vars)
	N = length(splits)
	H = splits[end][end]
	X = policy_rollout_only(env, x0, splits, vars)
	plot_rollout(X, Xref, splits)
	con_x, con_u = constraints(X, obj, splits, vars)
	for i = 1:N-1
		obj.λx[i] += obj.ρx * con_x[i]
	end
	for i = 1:H-1
		obj.λu[i] += obj.ρu * con_u[i]
	end
	return nothing
end

function linesearch(env, obj, Xref, splits, vars, Δ, f, vio)
	N = length(splits)
	α = 1.0
	for i = 1:10
		f1 = cost_only(env, obj, Xref, splits, vars+α*Δ)
		X1 = policy_rollout_only(env, x0, splits, vars+α*Δ)
		(f1 <= f || violation(X1, obj, splits, vars) <= vio) && break
		α *= 0.5
	end
	return α
end
