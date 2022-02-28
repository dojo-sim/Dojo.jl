"""
    SolverOptions{T}

    Options of the primal-dual interior point solver used to simulate the dynamics.

    rtol - defaults to 1e-6 - tolerance on residual violations (equality constraints)
    btol - defaults to 1e-3 - tolerance on bilinear violations (complementarity constraints)
    ls_scale - defaults to 0.5 - line search scaling factor (α_new ← ls_scale * α_current)
    max_iter - defaults to 50 - maximum number of Newton iterations
    max_ls - defaults to 10 - maximum number of line search steps
    undercut - defaults to Inf - complementarity slackness target; solver will aim at reaching complementarity violation = btol / undercut
    no_progress_max - defaults to 3 - number of Newton's iterations without progress trigerring the rescaling of the target complementarity violation
    no_progress_undercut - defaults to 10 - undercut scaling factor (target_new ← target_current / no_progress_undercut)
    verbose - defaults to false - printing the status of the solver during the solve
"""
@with_kw mutable struct SolverOptions{T}
    rtol::T=1.0e-6
    btol::T=1.0e-3
    ls_scale::T=0.5
    max_iter::Int=50
    max_ls=10
    undercut::T=Inf
    no_progress_max::Int=3
    no_progress_undercut::T=10.0
    verbose::Bool=false
end
