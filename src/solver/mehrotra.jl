# solver options
@with_kw struct InteriorPointOptions{T}
    rtol::T=1.0e-6   # residual violation tolerance (equality constraints)
    btol::T=1.0e-3   # bilinear violation tolerance (complementarity constraints)
    ls_scale::T=0.5  # line search scaling factor (α_new ← ls_scale * α_current)
    max_iter::Int=50 # maximum number of Newton iterations
    max_ls::Int=10   # maximum number of line search steps
    undercut::T=Inf  # complementarity slackness target; solver will aim at reaching κ_vio = btol / undercut
    verbose::Bool=false
end


# solver
function mehrotra!(mechanism::Mechanism; opts=InteriorPointOptions())
    system = mechanism.system
    eqcs = mechanism.eqconstraints
    bodies = mechanism.bodies
    ineqcs = mechanism.ineqconstraints

	resetVars!.(ineqcs, scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option
	resetVars!.(eqcs, scale=1.0) # resets the values of s and γ to the scaled neutral vector; TODO: solver option
    mechanism.μ = 0.0
	μtarget = 0.0
	no_progress = 0
    undercut = opts.undercut
    α = 1.0

	initial_state!.(ineqcs.values) # TODO: redundant with resetVars--remove
    setentries!(mechanism) # compute the residual

    bvio = bilinear_violation(mechanism) # does not require to apply setentries!
    rvio = residual_violation(mechanism) # does not require to apply setentries!

	opts.verbose && println("-----------------------------------------------------------------")

    for n = Base.OneTo(opts.max_iter)

        opts.verbose && solver_status(mechanism, α, rvio, bvio, n, μtarget, undercut)

        ((rvio < opts.rtol) && (bvio < opts.btol)) && break
		(n == opts.max_iter) && (opts.verbose && (@warn "failed mehrotra"))

        # affine search direction
		μ = 0.0
		pullresidual!(mechanism)                # store the residual inside mechanism.residual_entries
        ldu_factorization!(mechanism.system)    # factorize system, modifies the matrix in place
        pullmatrix!(mechanism)                  # store the factorized matrix inside mechanism.matrix_entries
        ldu_backsubstitution!(mechanism.system) # solve system, modifies the vector in place

		αaff = feasibilityStepLength!(mechanism; τort=0.95, τsoc=0.95, scaling=false) # uses system.vector_entries which holds the search drection
		ν, νaff = centering!(mechanism, αaff)
		σcentering = clamp(νaff / (ν + 1e-20), 0.0, 1.0)^3

		# corrected search direction
		μtarget = max(σcentering * ν, opts.btol / undercut)
		mechanism.μ = μtarget
		correction!(mechanism) # update the residual in mechanism.residual_entries
		mechanism.μ = 0.0

		pushresidual!(mechanism)                # cache residual + correction
        pushmatrix!(mechanism)                  # restore the factorized matrix
        ldu_backsubstitution!(mechanism.system) # solve system

		τ = max(0.95, 1 - max(rvio, bvio)^2) # τ = 0.95
		α = feasibilityStepLength!(mechanism; τort=τ, τsoc=min(τ, 0.95), scaling=false) # uses system.vector_entries which holds the corrected search drection

		# steps taken without making progress
		rvio_, bvio_ = lineSearch!(mechanism, α, rvio, bvio, opts; warning=false)
		made_progress = (!(rvio_ < opts.rtol) && (rvio_ < 0.8rvio)) || (!(bvio_ < opts.btol) && (bvio_ < 0.8bvio)) # we only care when progress is made while the tolerance is not met
		made_progress ? no_progress = max(no_progress - 1, 0) : no_progress += 1
		rvio, bvio = rvio_, bvio_
		(no_progress >= 3) && (undercut *= 10.0)

		# update solution
        foreach(updatesolution!, bodies)
        foreach(updatesolution!, eqcs)
        foreach(updatesolution!, ineqcs)

		# recompute Jacobian and residual
        setentries!(mechanism)
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
        "   bvio", scn(bvio),
        "   rvio", scn(rvio),
        "   α", scn(α),
        "   μ", scn(μtarget),
        "   |res|∞", scn(res),
        "   |Δ|∞", scn(Δvar),
        "   ucut", scn(undercut))
end

function initial_state!(ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    initial_state_ort(ineqc.γsol[1], ineqc.ssol[1])
    initial_state_ort(ineqc.γsol[2], ineqc.ssol[2])
    return nothing
end

function initial_state!(ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}}}
	γort, sort = initial_state_ort(ineqc.γsol[1][1:1], ineqc.ssol[1][1:1])
	γsoc, ssoc = initial_state_soc(ineqc.γsol[1][2:4], ineqc.ssol[1][2:4])
	ineqc.γsol[1] = [γort; γsoc]
	ineqc.ssol[1] = [sort; ssoc]
	γort, sort = initial_state_ort(ineqc.γsol[2][1:1], ineqc.ssol[2][1:1])
	γsoc, ssoc = initial_state_soc(ineqc.γsol[2][2:4], ineqc.ssol[2][2:4])
	ineqc.γsol[2] = [γort; γsoc]
	ineqc.ssol[2] = [sort; ssoc]
    return nothing
end

function initial_state_ort(γ::AbstractVector{T}, s::AbstractVector{T}; ϵ::T = 1e-20) where {T}
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

function initial_state_soc(γ::AbstractVector{T}, s::AbstractVector{T}; ϵ::T = 1e-20) where {T}
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
        component = getcomponent(mechanism, id)
        correction!(mechanism, residual_entries[id], getentry(system, id), component)
    end
	return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry,
		step_entry::Entry, component::Component)
    return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry,
		ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	μ = mechanism.μ
	residual_entry.value += [- Δs .* Δγ .+ μ; szeros(N½)]
    return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry,
		eqc::EqualityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
	cor = []
	for (i, joint) in enumerate(eqc.constraints)
        Δs, Δγ = get_sγ(joint, step_entry.value[λindex(eqc,i)])
		μ = mechanism.μ
		push!(cor, [- Δs .* Δγ .+ μ; szeros(length(joint) + 2 * joint_limits_length(joint))])
    end
	cor = vcat(cor...)
	residual_entry.value += cor
    return
end

@inline function correction!(mechanism::Mechanism, residual_entry::Entry, step_entry::Entry,
		ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}
	cont = ineqc.constraints[1]
	μ = mechanism.μ
	Δs = step_entry.value[1:N½]
    Δγ = step_entry.value[N½ .+ (1:N½)]
	residual_entry.value += [[-Δs[1] * Δγ[1]; -cone_product(Δs[2:4], Δγ[2:4])] + μ * neutral_vector(cont); szeros(N½)]
    return
end

function pullresidual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.residual_entries[i].value = mechanism.system.vector_entries[i].value
	end
	return
end

function pushresidual!(mechanism::Mechanism)
	for i in eachindex(mechanism.residual_entries)
		mechanism.system.vector_entries[i].value = mechanism.residual_entries[i].value
	end
	return
end

function pullmatrix!(mechanism::Mechanism)
	mechanism.matrix_entries.nzval .= mechanism.system.matrix_entries.nzval #TODO: make allocation free
	return
end

function pushmatrix!(mechanism::Mechanism)
	mechanism.system.matrix_entries.nzval .= mechanism.matrix_entries.nzval #TODO: make allocation free
	return
end

function savediagonalinverses!(mechanism::Mechanism)
	for i in eachindex(mechanism.diagonal_inverses)
		mechanism.diagonal_inverses[i] = deepcopy(mechanism.system.diagonal_inverses[i])
	end
	return
end

function pushdiagonalinverses!(mechanism::Mechanism)
	for i in eachindex(mechanism.diagonal_inverses)
		mechanism.system.diagonal_inverses[i] = deepcopy(mechanism.diagonal_inverses[i])
	end
	return
end
