function line_search!(mechanism::Mechanism, α, rvio, bvio, opts)
    scale = 0
    system = mechanism.system

    rvio_cand, bvio_cand = Inf * ones(2)
    for n = Base.OneTo(opts.max_ls)
        for contact in mechanism.contacts
            # candidate_step!(α[contact.id,:], contact, get_entry(system, contact.id), scale)
            # candidate_step!(fill(minimum(α[contact.id,:]), 2), contact, get_entry(system, contact.id), scale)
            candidate_step!([minimum(α[:,1]), minimum(α[:,2])] , contact, get_entry(system, contact.id), scale)
            # candidate_step!([minimum(α), minimum(α)] , contact, get_entry(system, contact.id), scale)
        end
        for joint in mechanism.joints
            # candidate_step!(α[joint.id,:], joint, get_entry(system, joint.id), scale)
            # candidate_step!(fill(minimum(α[joint.id,:]), 2), joint, get_entry(system, joint.id), scale)
            candidate_step!([minimum(α[:,1]), minimum(α[:,2])] , joint, get_entry(system, joint.id), scale)
            # candidate_step!([minimum(α), minimum(α)] , joint, get_entry(system, joint.id), scale)
        end
        for body in mechanism.bodies
            ϕmax = 3.9 / mechanism.timestep^2
            candidate_step!(mean(α), mechanism, body, get_entry(system, body.id), scale, ϕmax = ϕmax) # TODO this α is always one
            # candidate_step!(minimum(α), mechanism, body, get_entry(system, body.id), scale, ϕmax = ϕmax) # TODO this α is always one
            if dot(body.state.ϕsol[2], body.state.ϕsol[2]) > 3.91 / mechanism.timestep^2
                error("Excessive angular velocity. Body-ID: $(string(body.name)) " * string(body.id) * ", ω: " * string(body.state.ϕsol[2]) * ".")
            end
        end

        rvio_cand = residual_violation(mechanism)
        bvio_cand = bilinear_violation(mechanism)

        if (rvio_cand > rvio) && (bvio_cand > bvio)
            scale += 1
        else
            return rvio_cand, bvio_cand
        end
    end

    return rvio_cand, bvio_cand
end

function cone_line_search!(mechanism::Mechanism{T,Nn};
    τort::T=0.95,
    τsoc::T=0.95,
    scaling::Bool=false) where{T,Nn}

    system = mechanism.system

    α = ones(T,Nn,2)
    for contact in mechanism.contacts
        α[contact.id,:] = cone_line_search!(α[contact.id,:], mechanism, contact, get_entry(system, contact.id), τort, τsoc; scaling = scaling)
    end
    for joint in mechanism.joints
        α[joint.id,:] = cone_line_search!(α[joint.id,:], mechanism, joint, get_entry(system, joint.id), τort, τsoc; scaling = scaling)
    end

    return α
end

function cone_line_search!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc;
        scaling::Bool=false) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}

    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    αs_ort = positive_orthant_step_length(s[1:1], Δs[1:1], τ = τort)
    αγ_ort = positive_orthant_step_length(γ[1:1], Δγ[1:1], τ = τort)
    αs_soc = second_order_cone_step_length(s[2:4], Δs[2:4]; τ = τsoc)
    αγ_soc = second_order_cone_step_length(γ[2:4], Δγ[2:4]; τ = τsoc)

    α = [min(αs_ort, αs_soc), min(αγ_ort, αγ_soc)]
    return α
    # return min(α, αs_soc, αγ_soc, αs_ort, αγ_ort)
end

function cone_line_search!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc;
        scaling::Bool=false) where {T,N,Nc,Cs<:Union{ImpactContact{T,N},LinearContact{T,N}},N½}

    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]


    αs_ort = positive_orthant_step_length(s, Δs, τ = τort)
    αγ_ort = positive_orthant_step_length(γ, Δγ, τ = τort)

    return [αs_ort, αγ_ort]
    return min(α, αs_ort, αγ_ort)
end

function cone_line_search!(α, mechanism, joint::JointConstraint{T,N,Nc},
        vector_entry::Entry, τort, τsoc;
        scaling::Bool=false) where {T,N,Nc}

    αs_ort = 1.0
    αγ_ort = 1.0
    for (i, element) in enumerate([joint.translational, joint.rotational])
        s, γ = split_impulses(element, joint.impulses[2][joint_impulse_index(joint,i)])
        Δs, Δγ = split_impulses(element,  vector_entry.value[joint_impulse_index(joint,i)])

        αs_ort = min(αs_ort, positive_orthant_step_length(s, Δs, τ = τort))
        αγ_ort = min(αγ_ort, positive_orthant_step_length(γ, Δγ, τ = τort))
        # α = min(α, αs_ort, αγ_ort)
    end

    return [αs_ort, αγ_ort]
    # return α
end

function positive_orthant_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
    τ::T = 0.99) where T

    α = 1.0
    for i in eachindex(λ)
        if Δ[i] < 0 # safer
            α = min(α, - τ * λ[i] / Δ[i])
        end
    end

    return α
end

function second_order_cone_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T=0.99,
        ϵ::T=1e-14) where T

    # check Section 8.2 CVXOPT
    λ0 = λ[1]
    λ_λ = max(λ0^2 - λ[2:end]' * λ[2:end], 1e-25)
    if λ_λ < 0.0
        # @show λ_λ
        @warn "should always be positive"
    end
    λ_λ += ϵ
    λ_Δ = λ0 * Δ[1] - λ[2:end]' * Δ[2:end] + ϵ

    ρs = λ_Δ / λ_λ
    ρv = Δ[2:end] / sqrt(λ_λ)
    ρv -= (λ_Δ / sqrt(λ_λ) + Δ[1]) / (λ0 / sqrt(λ_λ) + 1) * λ[2:end] / λ_λ
    α = 1.0
    if norm(ρv) - ρs > 0.0
        α = min(α, τ / (norm(ρv) - ρs))
    end

    return α
end

function candidate_step!(α, mechanism::Mechanism, body::Body, vector_entry::Entry, scale; ϕmax = Inf)
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[1; 2; 3]] # TODO this α is always one
    body.state.ϕsol[2] = body.state.ϕsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[4; 5; 6]] # TODO this α is always one
    ϕ = body.state.ϕsol[2]
    ϕdot = dot(ϕ, ϕ)
    if ϕdot > ϕmax
        println("clipping ", scale, scn((ϕdot - ϕmax) / ϕmax), " ", scn(ϕdot), " ", scn(ϕmax), " ", body.name)
        body.state.ϕsol[2] *= ϕmax / ϕdot # this is overkill, but works better than sqrt(ϕmax/ϕdot)
    end
    return
end

function candidate_step!(α, joint::JointConstraint, vector_entry::Entry, scale)
    joint.impulses[2] = joint.impulses[1] + 1.0 / (2^scale) * minimum(α) * vector_entry.value # TODO maybe we need to split between s and γ
    return
end

function candidate_step!(α, contact::ContactConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, scale) where {T,N,Nc,Cs,N½}
    contact.impulses_dual[2] = contact.impulses_dual[1] + 1 / (2^scale) * α[1] * vector_entry.value[SVector{N½,Int64}(1:N½)]
    contact.impulses[2] = contact.impulses[1] + 1 / (2^scale) * α[2] * vector_entry.value[SVector{N½,Int64}(N½+1:N)]
    return
end
