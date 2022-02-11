@with_kw mutable struct SolverOptions{T} # TODO maybe remove the mutable, I used it for curriculum iLQR
    rtol::T=1.0e-6   # residual violation tolerance (equality constraints)
    btol::T=1.0e-3   # bilinear violation tolerance (complementarity constraints)
    ls_scale::T=0.5  # line search scaling factor (α_new ← ls_scale * α_current)
    max_iter::Int=50 # maximum number of Newton iterations
    max_ls::Int=10   # maximum number of line search steps
    undercut::T=Inf  # complementarity slackness target; solver will aim at reaching κ_vio = btol / undercut
    no_progress_max::Int=3 # number of iterations of no progress before rescaling undercut
    no_progress_undercut::T=10.0 # scaling for undercut if no progress is made
    verbose::Bool=false
end

function mehrotra!(mechanism::Mechanism; opts=SolverOptions())
	reset!.(mechanism.contacts, scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option
	reset!.(mechanism.joints, scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option

    mechanism.μ = 0.0
	μtarget = 0.0
	no_progress = 0
    undercut = opts.undercut
    α = 1.0

	initial_state!.(mechanism.contacts) # TODO: redundant with resetVars--remove
    set_entries!(mechanism) # compute the residual

    bvio = bilinear_violation(mechanism) # does not require to apply set_entries!
    rvio = residual_violation(mechanism) # does not require to apply set_entries!

	opts.verbose && println("-----------------------------------------------------------------")

    for n = Base.OneTo(opts.max_iter)

        opts.verbose && solver_status(mechanism, α, rvio, bvio, n, μtarget, undercut)

        ((rvio < opts.rtol) && (bvio < opts.btol)) && break
		(n == opts.max_iter) && (opts.verbose && (@warn "failed mehrotra"))

        # affine search direction
		μ = 0.0
		pull_residual!(mechanism)                # store the residual inside mechanism.residual_entries
        ldu_factorization!(mechanism.system)    # factorize system, modifies the matrix in place
        pull_matrix!(mechanism)                  # store the factorized matrix inside mechanism.matrix_entries
        ldu_backsubstitution!(mechanism.system) # solve system, modifies the vector in place

		αaff = feasibility_linesearch!(mechanism; τort=0.95, τsoc=0.95, scaling=false) # uses system.vector_entries which holds the search drection
		ν, νaff = centering!(mechanism, αaff)
		σcentering = clamp(νaff / (ν + 1e-20), 0.0, 1.0)^3

		# corrected search direction
		μtarget = max(σcentering * ν, opts.btol / undercut)
		mechanism.μ = μtarget
		correction!(mechanism) # update the residual in mechanism.residual_entries

		push_residual!(mechanism)                # cache residual + correction
        push_matrix!(mechanism)                  # restore the factorized matrix
        ldu_backsubstitution!(mechanism.system) # solve system

		τ = max(0.95, 1 - max(rvio, bvio)^2) # τ = 0.95
		α = feasibility_linesearch!(mechanism; τort=τ, τsoc=min(τ, 0.95), scaling=false) # uses system.vector_entries which holds the corrected search drection

		# steps taken without making progress
		rvio_, bvio_ = line_search!(mechanism, α, rvio, bvio, opts; warning=false)
		made_progress = (!(rvio_ < opts.rtol) && (rvio_ < 0.8rvio)) || (!(bvio_ < opts.btol) && (bvio_ < 0.8bvio)) # we only care when progress is made while the tolerance is not met
		made_progress ? no_progress = max(no_progress - 1, 0) : no_progress += 1
		rvio, bvio = rvio_, bvio_
		(no_progress >= opts.no_progress_max) && (undercut *= opts.no_progress_undercut)

		# update solution
        update_solution!.(mechanism.bodies)
        update_solution!.(mechanism.joints)
        update_solution!.(mechanism.contacts)

		# recompute Jacobian and residual
        set_entries!(mechanism)
    end

    return
end

function solver_status(mechanism::Mechanism, α, rvio, bvio, n, μtarget, undercut)
    fv = full_vector(mechanism.system)
    Δvar = norm(fv, Inf)
    fM = full_matrix(mechanism.system)
    fΔ = fM \ fv
    Δalt = norm(fΔ, Inf)
    ##################
    res = norm(fv, Inf)
    println(
        "n ", n,
        "   bvio", scn(bvio, digits=0),
        "   rvio", scn(rvio, digits=0),
        "   α", scn(α, digits=0),
        "   μ", scn(μtarget, digits=0),
        "   |res|∞", scn(res, digits=0),
        "   |Δ|∞", scn(Δvar, digits=0),
        # "   ucut", scn(undercut),
        )
end

function initial_state!(contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    initialize_positive_orthant!(contact.impulses[1], contact.impulses_dual[1])
    initialize_positive_orthant!(contact.impulses[2], contact.impulses_dual[2])
    return nothing
end

function initial_state!(contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:NonlinearContact{T,N}}
	γort, sort = initialize_positive_orthant!(contact.impulses[1][1:1], contact.impulses_dual[1][1:1])
	γsoc, ssoc = initialize_second_order_cone!(contact.impulses[1][2:4], contact.impulses_dual[1][2:4])
	contact.impulses[1] = [γort; γsoc]
	contact.impulses_dual[1] = [sort; ssoc]
	γort, sort = initialize_positive_orthant!(contact.impulses[2][1:1], contact.impulses_dual[2][1:1])
	γsoc, ssoc = initialize_second_order_cone!(contact.impulses[2][2:4], contact.impulses_dual[2][2:4])
	contact.impulses[2] = [γort; γsoc]
	contact.impulses_dual[2] = [sort; ssoc]
    return nothing
end

function initialize_positive_orthant!(γ::AbstractVector{T}, s::AbstractVector{T}; ϵ::T = 1e-20) where T
    δs = max(-1.5 * minimum(s), 0)
    δγ = max(-1.5 * minimum(γ), 0)

    sh = s .+ δs
    γh = γ .+ δγ
    δhs = 0.5 * transpose(sh) * γh / (sum(γh) + ϵ)
    δhγ = 0.5 * transpose(sh) * γh / (sum(sh) + ϵ)
    s0 = sh .+ δhs
    γ0 = γh .+ δhγ
	return γ0, s0
end

function initialize_second_order_cone!(γ::AbstractVector{T}, s::AbstractVector{T}; ϵ::T = 1e-20) where T
    e = [1.0; zeros(length(γ) - 1)] # identity element
    δs = max(-1.5 * (s[1] - norm(s[2:end])), 0)
    δγ = max(-1.5 * (γ[1] - norm(γ[2:end])), 0)

    sh = s + δs * e
    γh = γ + δγ * e
    δhs = 0.5 * transpose(sh) * γh / ((γh[1] + norm(γh[2,end])) + ϵ)
    δhγ = 0.5 * transpose(sh) * γh / ((sh[1] + norm(sh[2,end])) + ϵ)

    s0 = sh + δhs * e
    γ0 = γh + δhγ * e
	return γ0, s0
end

function correction!(mechanism)
	system = mechanism.system
	residual_entries = mechanism.residual_entries

    for id in reverse(system.dfs_list)
        node = get_node(mechanism, id)
        correction!(mechanism, residual_entries[id], get_entry(system, id), node)
    end
	return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, node::Node)
    return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, ::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	μ = mechanism.μ
	residual_entry.value += [- Δs .* Δγ .+ μ; szeros(N½)]
    return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}
	cont = contact.model
	μ = mechanism.μ
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	residual_entry.value += [[-Δs[1] * Δγ[1]; -cone_product(Δs[2:4], Δγ[2:4])] + μ * neutral_vector(cont); szeros(N½)]
    return
end

@inline function correction!(mechanism::Mechanism{T}, residual_entry::Entry, step_entry::Entry, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
	cor = correction(mechanism, step_entry, joint)
	residual_entry.value += cor
    return
end

@generated function correction(mechanism::Mechanism{T}, step_entry::Entry, joint::JointConstraint{T,N,Nc}) where {T,N,Nc}
    cor = [:(correction([joint.translational, joint.rotational][$i], step_entry.value[λindex(joint, $i)], mechanism.μ)) for i = 1:Nc]
    return :(vcat($(cor...)))
end

@inline function correction(joint::Joint{T,Nλ,Nb,N}, Δ, μ) where {T,Nλ,Nb,N}
    Δs, Δγ = get_sγ(joint, Δ)
	return [- Δs .* Δγ .+ μ; szeros(Nb + Nλ)]
end

function pull_residual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.residual_entries[i].value = mechanism.system.vector_entries[i].value
	end
	return
end

function push_residual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.system.vector_entries[i].value = mechanism.residual_entries[i].value
	end
	return
end

function pull_matrix!(mechanism::Mechanism)
	mechanism.matrix_entries.nzval .= mechanism.system.matrix_entries.nzval #TODO: make allocation free
	return
end

function push_matrix!(mechanism::Mechanism)
	mechanism.system.matrix_entries.nzval .= mechanism.matrix_entries.nzval #TODO: make allocation free
	return
end

function save_diagonal_inverses!(mechanism::Mechanism)
	for i in eachindex(mechanism.diagonal_inverses)
		mechanism.diagonal_inverses[i] = deepcopy(mechanism.system.diagonal_inverses[i])
	end
	return
end

function push_diagonal_inverses!(mechanism::Mechanism)
	for i in eachindex(mechanism.diagonal_inverses)
		mechanism.system.diagonal_inverses[i] = deepcopy(mechanism.diagonal_inverses[i])
	end
	return
end
