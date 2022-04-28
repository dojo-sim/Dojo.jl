"""
    augmented Lagrangian solve
"""
function IterativeLQR.constrained_ilqr_solve!(solver::Solver; callback!::Function=x->nothing)

	# verbose && printstyled("Iterative LQR\n",
	# 	color=:red, bold=true)

    # reset solver cache
    IterativeLQR.reset!(solver.data)

    # reset duals
    for (t, λ) in enumerate(solver.problem.objective.costs.constraint_dual)
        fill!(λ, 0.0)
	end

	# initialize penalty
	for (t, ρ) in enumerate(solver.problem.objective.costs.constraint_penalty)
        fill!(ρ, solver.options.initial_constraint_penalty)
	end

	for i = 1:solver.options.max_dual_updates
		solver.options.verbose && println("  al iter: $i")

		# primal minimization
		IterativeLQR.ilqr_solve!(solver)

		# update trajectories
		IterativeLQR.cost!(solver.data, solver.problem,
            mode=:nominal)

        # constraint violation
		solver.data.max_violation[1] <= solver.options.constraint_tolerance && break

        # dual ascent
		IterativeLQR.augmented_lagrangian_update!(solver.problem.objective.costs,
			scaling_penalty=solver.options.scaling_penalty,
            max_penalty=solver.options.max_penalty)

		# user-defined callback (continuation methods on the models etc.)
		callback!(solver)
	end

    return nothing
end
