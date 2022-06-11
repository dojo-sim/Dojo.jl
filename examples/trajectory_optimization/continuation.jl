function reset!(dynamics_model::Simulator; r_tol=1e-5, κ_tol = 1e-3, undercut = 5.0)
    dynamics_model.ip.opts.r_tol = r_tol
    dynamics_model.ip.opts.κ_tol = κ_tol
    dynamics_model.ip.opts.undercut = undercut
    return nothing
end

function continuation_callback!(solver::Solver, dymamics_model::RoboDojo.Simulator)
    dynamics_model.ip.opts.r_tol = max(dynamics_model.ip.opts.r_tol/3, 1e-7)
    dynamics_model.ip.opts.κ_tol = max(dynamics_model.ip.opts.κ_tol/3, 1e-5)
    println("r_tol $(scn(dynamics_model.ip.opts.r_tol))  " *
        "κ_tol $(scn(dynamics_model.ip.opts.κ_tol))")

    # check that we are not too far from a solution
    solver_kickout = solver.data.max_violation[1] > 1e3
    return solver_kickout
end
