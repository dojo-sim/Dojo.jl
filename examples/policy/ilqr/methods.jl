function continuation_callback!(solver::Solver, dymamics_model::RoboDojo.Simulator)
    dynamics_model.ip.opts.r_tol = max(dynamics_model.ip.opts.r_tol/3, 1e-7)
    dynamics_model.ip.opts.κ_tol = max(dynamics_model.ip.opts.κ_tol/3, 1e-5)
    println("r_tol $(scn(dynamics_model.ip.opts.r_tol))  " *
        "κ_tol $(scn(dynamics_model.ip.opts.κ_tol))")
    return nothing
end
