"""
    SolverOptions{T} 

    options and tolerances for interior-point solver 

    rtol: residual violation tolerance (equality constraints)
    btol: bilinear violation tolerance (complementarity constraints)
    ls_scale: line search scaling factor (α_new ← ls_scale * α_current)
    max_iter: maximum number of Newton iterations
    max_ls: maximum number of line search steps
    undercut: complementarity slackness target; solver will aim at reaching κ_vio = btol / undercut
    no_progress_max: number of iterations of no progress before rescaling undercut
    no_progress_undercut: scaling for undercut if no progress is made
    verbose: flag for printing
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