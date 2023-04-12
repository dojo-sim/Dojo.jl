function line_search!(mechanism::Mechanism, α, rvio, bvio, opts)
    scale = 0
    system = mechanism.system

    rvio_cand = Inf
    bvio_cand = Inf
    for n = Base.OneTo(opts.max_ls)
        for contact in mechanism.contacts
            candidate_step!(α, contact, get_entry(system, contact.id), scale)
        end
        for joint in mechanism.joints
            candidate_step!(α, joint, get_entry(system, joint.id), scale)
        end
        for body in mechanism.bodies
            ωmax = 3.9 / mechanism.timestep^2
            candidate_step!(α, mechanism, body, get_entry(system, body.id), scale, 
                ωmax=ωmax, verbose=opts.verbose)
            if dot(body.state.ωsol[2], body.state.ωsol[2]) > 3.91 / mechanism.timestep^2
                error("Excessive angular velocity. Body-ID: $(string(body.name)) " * string(body.id) * ", ω: " * string(body.state.ωsol[2]) * ".")
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

function cone_line_search!(mechanism::Mechanism;
    τort::T=0.95,
    τsoc::T=0.95) where T

    system = mechanism.system

    α = 1.0
    for contact in mechanism.contacts
        α = cone_line_search!(α, mechanism, contact, get_entry(system, contact.id), τort, τsoc)
    end
    for joint in mechanism.joints
        α = cone_line_search!(α, mechanism, joint, get_entry(system, joint.id), τort, τsoc)
    end

    return α
end

function cone_line_search!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc) where {T,N,Nc,Cs<:NonlinearContact{T,N},N½}

    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    αs_ort = positive_orthant_step_length(s[1:1], Δs[1:1], τ = τort)
    αγ_ort = positive_orthant_step_length(γ[1:1], Δγ[1:1], τ = τort)
    αs_soc = second_order_cone_step_length(s[2:4], Δs[2:4]; τ = τsoc)
    αγ_soc = second_order_cone_step_length(γ[2:4], Δγ[2:4]; τ = τsoc)

    return min(α, αs_soc, αγ_soc, αs_ort, αγ_ort)
end

function cone_line_search!(α, mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc) where {T,N,Nc,Cs<:Union{ImpactContact{T,N},LinearContact{T,N}},N½}

    s = contact.impulses_dual[2]
    γ = contact.impulses[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]


    αs_ort = positive_orthant_step_length(s, Δs, τ = τort)
    αγ_ort = positive_orthant_step_length(γ, Δγ, τ = τort)

    return min(α, αs_ort, αγ_ort)
end

function cone_line_search!(α, mechanism, joint::JointConstraint{T,N,Nc},
        vector_entry::Entry, τort, τsoc) where {T,N,Nc}

    
    s, γ = split_impulses(joint.translational, joint.impulses[2][joint_impulse_index(joint,1)])
    Δs, Δγ = split_impulses(joint.translational,  vector_entry.value[joint_impulse_index(joint,1)])
    αs_ort = positive_orthant_step_length(s, Δs, τ = τort)
    αγ_ort = positive_orthant_step_length(γ, Δγ, τ = τort)
    α = min(α, αs_ort, αγ_ort)

    s, γ = split_impulses(joint.rotational, joint.impulses[2][joint_impulse_index(joint,2)])
    Δs, Δγ = split_impulses(joint.rotational,  vector_entry.value[joint_impulse_index(joint,2)])
    αs_ort = positive_orthant_step_length(s, Δs, τ = τort)
    αγ_ort = positive_orthant_step_length(γ, Δγ, τ = τort)
    α = min(α, αs_ort, αγ_ort)

    return α
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

function candidate_step!(α, mechanism::Mechanism, body::Body, vector_entry::Entry, scale; 
    ωmax = Inf, verbose=true)
    body.state.vsol[2] = body.state.vsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[1; 2; 3]]
    body.state.ωsol[2] = body.state.ωsol[1] + 1 / (2^scale) * α * vector_entry.value[SA[4; 5; 6]]
    ω = body.state.ωsol[2]
    ωdot = dot(ω, ω)
    if ωdot > ωmax
        verbose && println("clipping ", scale, scn((ωdot - ωmax) / ωmax), " ", scn(ωdot), " ", scn(ωmax), " ", body.name)
        body.state.ωsol[2] *= ωmax / ωdot # this is overkill, but works better than sqrt(ωmax/ωdot)
    end
    return
end

function candidate_step!(α, joint::JointConstraint, vector_entry::Entry, scale)
    joint.impulses[2] = joint.impulses[1] + 1.0 / (2^scale) * α * vector_entry.value
    return
end

function candidate_step!(α::T, contact::ContactConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, scale) where {T,N,Nc,Cs,N½}
    contact.impulses_dual[2] = contact.impulses_dual[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(1:N½)]
    contact.impulses[2] = contact.impulses[1] + 1 / (2^scale) * α * vector_entry.value[SVector{N½,Int64}(N½+1:N)]
    return
end
