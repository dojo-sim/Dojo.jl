@with_kw mutable struct SolverOptions{T}
    rtol::T=1.0e-6   # residual violation tolerance (equality constraints)
    btol::T=1.0e-3   # bilinear violation tolerance (complementarity constraints)
    ls_scale::T=0.5  # line search scaling factor (α_new ← ls_scale * α_current)
    max_iter::Int=50 # maximum number of Newton iterations
    max_ls=10   # maximum number of line search steps
    undercut::T=Inf  # complementarity slackness target; solver will aim at reaching κ_vio = btol / undercut
    no_progress_max::Int=3 # number of iterations of no progress before rescaling undercut
    no_progress_undercut::T=10.0 # scaling for undercut if no progress is made
    verbose::Bool=false
end