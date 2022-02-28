"""
    SolverOptions{T}

    Options and tolerances of primal-dual interior point solver.

    rtol: tolerance on residual violations (equality constraints); defaults to 1e-6
    btol: tolerance on bilinear violations (complementarity constraints); defaults to 1e-3
    ls_scale: line search scaling factor (α_new ← ls_scale * α_current); defaults to 0.5
    max_iter: maximum number of Newton iterations; defaults to 50
    max_ls: maximum number of line search steps; defaults to 10
    undercut: complementarity slackness target; solver will aim at reaching complementarity violation = btol / undercut; defaults to Inf
    no_progress_max: number of Newton's iterations without progress trigerring the rescaling of the target complementarity violation; defaults to 3
    no_progress_undercut: undercut scaling factor (target_new ← target_current / no_progress_undercut); defaults to 10
    verbose: flag for printing the status of the solver during the solve; defaults to false
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
