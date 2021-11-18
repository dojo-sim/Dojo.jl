@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, component::Component)
    matrix_entry.value = ∂g∂ʳself(mechanism, component)
    vector_entry.value = -g(mechanism, component)
    return
end

@inline function setLU!(mechanism::Mechanism, matrix_entry_L::Entry, matrix_entry_U::Entry, componenta::Component, componentb::Component)
    L, U = ∂gab∂ʳba(mechanism, componenta, componentb)
    matrix_entry_L.value = L
    matrix_entry_U.value = U
    return
end

@inline function zeroLU!(matrix_entry_L::Entry, matrix_entry_U::Entry)
    matrix_entry_L.value *= 0
    matrix_entry_U.value *= 0
    return
end


function feasibilityStepLength!(mechanism::Mechanism; τort::T = 0.95, τsoc::T = 0.95) where {T}
    system = mechanism.system

    mechanism.α = 1.0

    for ineqc in mechanism.ineqconstraints
        feasibilityStepLength!(mechanism, ineqc, getentry(system, ineqc.id), τort, τsoc)
    end

    return
end

function feasibilityStepLength!(mechanism::Mechanism{T}, body::Body{T}, vector_entry::Entry; τ::T = 0.99) where {T}
    Δt = mechanism.Δt
    Δω = vector_entry.value[SA[4; 5; 6]]
    ω = body.state.ωsol[2]

    α = polynomial_step_length(ω, Δω, τ * 1/Δt)

    (α < mechanism.α) && (mechanism.α = α)
    return
end

function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}

    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]

    αs_ort = 1.0
    αγ_ort = 1.0
    for i = 1:1
        if Δs[i] < 0 # safer
            αs_ort = min(αs_ort, - τort * s[i] / Δs[i])
        end
        if Δγ[i] < 0 # safer
            αγ_ort = min(αγ_ort, - τort * γ[i] / Δγ[i])
        end
    end

    αs = soc_step_length(s[2:4], Δs[2:4]; τ = τsoc)
    αγ = soc_step_length(γ[2:4], Δγ[2:4]; τ = τsoc)
    αγ < 1e-6 && println("γ: ", scn.(γ[2:4]))
    αγ < 1e-6 && println("Δγ:", scn.(Δγ[2:4]))
    α = min(αs, αγ, αs_ort, αγ_ort)
    (α > 0) && (α < mechanism.α) && (mechanism.α = α)
    return
end

function feasibilityStepLength!(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc) where {T,N,Nc,Cs<:Tuple{Union{ImpactBound11{T,N},LinearContactBound11{T,N}}},N½}
    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]

    αs_ort = 1.0
    αγ_ort = 1.0
    for i = 1:N½
        if Δs[i] < 0 # safer
            αs_ort = min(αs_ort, - τort * s[i] / Δs[i])
        end
        if Δγ[i] < 0 # safer
            αγ_ort = min(αγ_ort, - τort * γ[i] / Δγ[i])
        end
    end
    α = min(αs_ort, αγ_ort)
    (α > 0) && (α < mechanism.α) && (mechanism.α = α)
    return
end

function centering!(mechanism::Mechanism, αaff::T) where {T}
    system = mechanism.system

    n = 0
    mechanism.ν = 0.0
    mechanism.νaff = 0.0
    for ineqc in mechanism.ineqconstraints
        n += centering!(mechanism, ineqc, getentry(system, ineqc.id), αaff)
    end

    mechanism.ν /= n
    mechanism.νaff /= n
    return
end

function centering!(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, αaff) where {T,N,Nc,Cs,N½}
    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    mechanism.ν += dot(s, γ)
    mechanism.νaff += dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    return N½
end

function soc_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T = 0.99, ϵ::T = 1e-14) where {T}
    # check Section 8.2 CVXOPT
    # The CVXOPT linear and quadratic cone program solvers

    # Adding to slack ϵ to make sure that we never get out of the cone
    λ0 = λ[1]
    λ_λ = max(λ0^2 - λ[2:end]' * λ[2:end], 1e-25)
    if λ_λ < 0.0
        @show λ_λ
        @warn "should always be positive"
    end
    λ_λ += ϵ
    λ_Δ = λ0 * Δ[1] - λ[2:end]' * Δ[2:end] + ϵ

    ρs = λ_Δ / λ_λ
    ρv = Δ[2:end] / sqrt(λ_λ)
    ρv -= (λ_Δ / sqrt(λ_λ) + Δ[1]) / (λ0 / sqrt(λ_λ) + 1) * λ[2:end] / λ_λ
    # we make sre that the inverse always exists with ϵ,
    # if norm(ρv) - ρs) is negative (Δ is pushing towards a more positive cone)
        # the computation is ignored and we get the maximum value for α = 1.0
    # else we have α = τ / norm(ρv) - ρs)
    # we add ϵ to the denumerator to ensure strict positivity and avoid 1e-16 errors.
    α = 1.0
    if norm(ρv) - ρs > 0.0
        α = min(α, τ / (norm(ρv) - ρs))
    end
    return α
end

function setentries!(mechanism::Mechanism)
    system = mechanism.system

    for id in reverse(system.dfs_list)
        for childid in system.cyclic_children[id]
            zeroLU!(getentry(system, id, childid), getentry(system, childid, id))
        end

        component = getcomponent(mechanism, id)
        setDandΔs!(mechanism, getentry(system, id, id), getentry(system, id), component)

        for childid in children(system,id)
            setLU!(mechanism, getentry(system, id, childid), getentry(system, childid, id), component, getcomponent(mechanism, childid))
        end
    end

    return
end

@inline function updatesolution!(body::Body)
    body.state.vsol[1] = body.state.vsol[2]
    body.state.ωsol[1] = body.state.ωsol[2]
    return
end

@inline function updatesolution!(eqc::EqualityConstraint)
    eqc.λsol[1] = eqc.λsol[2]
    return
end

@inline function updatesolution!(ineqc::InequalityConstraint)
    ineqc.ssol[1] = ineqc.ssol[2]
    ineqc.γsol[1] = ineqc.γsol[2]
    return
end

@inline function normΔs(body::Body)
    d1 = body.state.vsol[2] - body.state.vsol[1]
    d2 = body.state.ωsol[2] - body.state.ωsol[1]
    return dot(d1, d1) + dot(d2, d2)
end

@inline function normΔs(eqc::EqualityConstraint)
    d = eqc.λsol[2] - eqc.λsol[1]
    return dot(d, d)
end

@inline function normΔs(ineqc::InequalityConstraint)
    d1 = ineqc.ssol[2] - ineqc.ssol[1]
    d2 = ineqc.γsol[2] - ineqc.γsol[1]
    return dot(d1, d1) + dot(d2, d2)
end

function f(x, a, b, c)
    return a * x^2  + b * x + c
end

function polynomialsolution(a, b, c)
    Δ = b^2 - 4a*c
    xp = (-b + sqrt(Δ))/ (2a)
    xm = (-b - sqrt(Δ))/ (2a)
    return xp, xm
end

function polynomial_step_length(ω::AbstractVector{T}, Δ::AbstractVector{T}, β::T) where {T}
    a = Δ' * Δ
    b = 2Δ' * ω
    c = ω' * ω - β^2
    α = maximum(polynomialsolution(a, b, c))
    α = clamp(α, 0.0, 1.0)
    return α
end
