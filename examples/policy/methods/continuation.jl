function reset!(env::Environment; rtol=1e-4, btol=1e-2, undercut=5.0)
    env.opts_step.rtol = rtol
    env.opts_step.btol = btol
    env.opts_step.undercut = undercut
    env.opts_grad.rtol = rtol
    env.opts_grad.btol = btol
    env.opts_grad.undercut = undercut
    return nothing
end

function continuation_callback!(solver::Solver, env::Environment; ρ=3, build::Bool=false)
    # contact smoothness continuation
    env.opts_step.rtol = max(1e-6, env.opts_step.rtol/ρ)
    env.opts_step.btol = max(1e-4, env.opts_step.btol/ρ)
    env.opts_grad.rtol = max(1e-6, env.opts_grad.rtol/ρ)
    env.opts_grad.btol = max(1e-4, env.opts_grad.btol/ρ)

    # visualize current policy
    u = solver.problem.actions
    x = IterativeLQR.rollout(model, x1, u)
    DojoEnvironments.visualize(env, x, build=build)

    println("r_tol $(scn(env.opts_grad.rtol))  " * "κ_tol $(scn(env.opts_grad.btol))")
    return nothing
end
