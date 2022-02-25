function mehrotra!(mechanism::Mechanism; opts=SolverOptions())
	reset!.(mechanism.contacts, scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option
	reset!.(mechanism.joints,   scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option

	status = :failed
    mechanism.μ = 0.0
	μtarget = 0.0
	no_progress = 0
    undercut = opts.undercut
    α = 1.0

	initialize!.(mechanism.contacts) # TODO: redundant with resetVars--remove
    set_entries!(mechanism) # compute the residual

    bvio = bilinear_violation(mechanism) # does not require to apply set_entries!
    rvio = residual_violation(mechanism) # does not require to apply set_entries!
	opts.verbose && solver_header()
    for n = Base.OneTo(opts.max_iter)
        opts.verbose && solver_status(mechanism, α, rvio, bvio, n, μtarget, undercut)

        ((rvio < opts.rtol) && (bvio < opts.btol)) && (status=:success; break)
		(n == opts.max_iter) && (opts.verbose && (@warn "failed mehrotra"))

        # affine search direction
		μ = 0.0
		pull_residual!(mechanism)               # store the residual inside mechanism.residual_entries
        ldu_factorization!(mechanism.system)    # factorize system, modifies the matrix in place
        pull_matrix!(mechanism)                 # store the factorized matrix inside mechanism.matrix_entries
        ldu_backsubstitution!(mechanism.system) # solve system, modifies the vector in place

		αaff = cone_line_search!(mechanism; τort=0.95, τsoc=0.95, scaling=false) # uses system.vector_entries which holds the search drection
		ν, νaff = centering!(mechanism, αaff)
		σcentering = clamp(νaff / (ν + 1e-20), 0.0, 1.0)^3

		# corrected search direction
		μtarget = max(σcentering * ν, opts.btol / undercut)
		mechanism.μ = μtarget
		correction!(mechanism) # update the residual in mechanism.residual_entries

		push_residual!(mechanism)               # cache residual + correction
        push_matrix!(mechanism)                 # restore the factorized matrix
        ldu_backsubstitution!(mechanism.system) # solve system

		τ = max(0.95, 1 - max(rvio, bvio)^2) # τ = 0.95
		α = cone_line_search!(mechanism; τort=τ, τsoc=min(τ, 0.95), scaling=false) # uses system.vector_entries which holds the corrected search drection

		# steps taken without making progress
		rvio_, bvio_ = line_search!(mechanism, α, rvio, bvio, opts)
		made_progress = (!(rvio_ < opts.rtol) && (rvio_ < 0.8rvio)) || (!(bvio_ < opts.btol) && (bvio_ < 0.8bvio)) # we only care when progress is made while the tolerance is not met
		made_progress ? no_progress = max(no_progress - 1, 0) : no_progress += 1
		rvio, bvio = rvio_, bvio_
		(no_progress >= opts.no_progress_max) && (undercut *= opts.no_progress_undercut)

		# update solution
        update!.(mechanism.bodies)
        update!.(mechanism.joints)
        update!.(mechanism.contacts)

		# recompute Jacobian and residual
        set_entries!(mechanism)
    end

    return status
end

function solver_header()
	println("                                                 ")
	println("n    bvio    rvio     α       μ     |res|∞   |Δ|∞")
	println("–––––––––––––––––––––––––––––––––––––––––––––––––")
end

function status(mechanism::Mechanism, α, rvio, bvio, n, μtarget, undercut)
    fv = full_vector(mechanism.system)
    Δvar = norm(fv, Inf)
    fM = full_matrix(mechanism.system)
    fΔ = fM \ fv
    Δalt = norm(fΔ, Inf)
    res = norm(fv, Inf)
	println(
        n,
        "   ", scn(bvio, digits=0),
        "   ", scn(rvio, digits=0),
        "   ", scn(α, digits=0),
        "   ", scn(μtarget, digits=0),
        "   ", scn(res, digits=0),
        "   ", scn(Δvar, digits=0),
        # "   ucut", scn(undercut),
        )
end

function initial_state!(contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    initialize_positive_orthant!(contact.impulses[1], contact.impulses_dual[1])
    initialize_positive_orthant!(contact.impulses[2], contact.impulses_dual[2])
    return nothing
end
