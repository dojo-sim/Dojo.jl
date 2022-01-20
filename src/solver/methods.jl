@inline function setDandΔs!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry, component::Component)
    matrix_entry.value = ∂g∂z(mechanism, component)
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

#######


## Complementarity
function complementarity(mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}; scaling::Bool = false) where {T,N,Nc,Cs}
    c = []
    for (i, joint) in enumerate(eqc.constraints)
        λi = eqc.λsol[2][λindex(eqc, i)]
        si, γi = get_sγ(joint, λi)
        push!(c, si .* γi)
    end
    return vcat(c...)
end


function feasibilityStepLength!(mechanism::Mechanism; τort::T=0.95, τsoc::T=0.95, scaling::Bool=false) where {T}
    system = mechanism.system

    α = 1.0
    for ineqc in mechanism.ineqconstraints
        α = feasibilityStepLength!(α, mechanism, ineqc, getentry(system, ineqc.id), τort, τsoc; scaling = scaling)
    end
    for eqc in mechanism.eqconstraints
        α = feasibilityStepLength!(α, mechanism, eqc, getentry(system, eqc.id), τort, τsoc; scaling = scaling)
    end
    return α
end

function feasibilityStepLength!(α, mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc; scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{ContactBound{T,N}},N½}

    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    αs_ort = ort_step_length(s[1:1], Δs[1:1], τ = τort)
    αγ_ort = ort_step_length(γ[1:1], Δγ[1:1], τ = τort)
    αs_soc = soc_step_length(s[2:4], Δs[2:4]; τ = τsoc)
    αγ_soc = soc_step_length(γ[2:4], Δγ[2:4]; τ = τsoc)

    return min(α, αs_soc, αγ_soc, αs_ort, αγ_ort)
end

function feasibilityStepLength!(α, mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½},
        vector_entry::Entry, τort, τsoc; scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{Union{ImpactBound{T,N},LinearContactBound{T,N}}},N½}
    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]


    αs_ort = ort_step_length(s, Δs, τ = τort)
    αγ_ort = ort_step_length(γ, Δγ, τ = τort)

    return min(α, αs_ort, αγ_ort)
end

function feasibilityStepLength!(α, mechanism, eqc::EqualityConstraint{T,N,Nc,Cs},
        vector_entry::Entry, τort, τsoc; scaling::Bool = false) where {T,N,Nc,Cs}

    for (i, joint) in enumerate(eqc.constraints)
        s, γ = get_sγ(joint, eqc.λsol[2][λindex(eqc,i)])
        Δs, Δγ = get_sγ(joint,  vector_entry.value[λindex(eqc,i)])

        αs_ort = ort_step_length(s, Δs, τ = τort)
        αγ_ort = ort_step_length(γ, Δγ, τ = τort)
        α = min(α, αs_ort, αγ_ort)
    end

    return α
end

function centering!(mechanism::Mechanism, αaff::T) where {T}
    system = mechanism.system

    n = 0
    ν = 0.0
    νaff = 0.0
    for ineqc in mechanism.ineqconstraints
        ν, νaff, n = centering!(ν, νaff, n, mechanism, ineqc, getentry(system, ineqc.id), αaff)
    end

    for eqc in mechanism.eqconstraints
        ν, νaff, n = centering!(ν, νaff, n, mechanism, eqc, getentry(system, eqc.id), αaff)
    end

    ν /= n
    νaff /= n
    return ν, νaff
end

function centering!(ν, νaff, n, mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}, vector_entry::Entry, αaff) where {T,N,Nc,Cs,N½}
    s = ineqc.ssol[2]
    γ = ineqc.γsol[2]
    Δs = vector_entry.value[1:N½]
    Δγ = vector_entry.value[N½ .+ (1:N½)]
    ν += dot(s, γ)
    νaff += dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
    n += cone_degree(ineqc)
    return ν, νaff, n
end

function centering!(ν, νaff, n, mechanism, eqc::EqualityConstraint{T,N,Nc,Cs}, vector_entry::Entry, αaff) where {T,N,Nc,Cs}
    for (i, joint) in enumerate(eqc.constraints)
        s, γ = get_sγ(joint, eqc.λsol[2][λindex(eqc,i)])
        Δs, Δγ = get_sγ(joint, vector_entry.value[λindex(eqc,i)])
        ν += dot(s, γ)
        νaff += dot(s + αaff * Δs, γ + αaff * Δγ) # plus or minus
        n += length(s)
    end
    return ν, νaff, n
end

function ort_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T}; τ::T = 0.99) where {T}
    α = 1.0
    for i in eachindex(λ)
        if Δ[i] < 0 # safer
            α = min(α, - τ * λ[i] / Δ[i])
        end
    end
    return α
end

function soc_step_length(λ::AbstractVector{T}, Δ::AbstractVector{T};
        τ::T=0.99, ϵ::T=1e-14) where {T}
    # check Section 8.2 CVXOPT
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
    body.state.ϕsol[1] = body.state.ϕsol[2]
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
    d2 = body.state.ϕsol[2] - body.state.ϕsol[1]
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

@inline function residual_violation(mechanism::Mechanism)
    violation = 0.0
    for eq in mechanism.eqconstraints
        res = g(mechanism, eq)
        shift = 0
        for (i, joint) in enumerate(eq.constraints)
            Nλ = λlength(joint)
            Nb = blength(joint)
            subres = res[shift + 2Nb .+ (1:Nλ)]
            violation = max(violation, norm(subres, Inf))
            shift += ηlength(joint)
        end
    end
    for body in mechanism.bodies
        res = g(mechanism, body)
        violation = max(violation, norm(res, Inf))
    end
    for ineq in mechanism.ineqconstraints
        res = g(mechanism, ineq)
        violation = max(violation, norm(res, Inf))
    end
    return violation
end

@inline function bilinear_violation(mechanism::Mechanism)
    violation = 0.0
    for ineqc in mechanism.ineqconstraints
        comp = complementarity(mechanism, ineqc)
        violation = max(violation, norm(comp, Inf))
    end
    for eqc in mechanism.eqconstraints
        comp = complementarity(mechanism, eqc)
        violation = max(violation, norm(comp, Inf))
    end
    return violation
end
