function initialize!(contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs}
    initialize_positive_orthant!(contact.impulses[1], contact.impulses_dual[1])
    initialize_positive_orthant!(contact.impulses[2], contact.impulses_dual[2])
    return nothing
end

function initialize!(contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:NonlinearContact{T,N}}
	γort, sort = initialize_positive_orthant!(contact.impulses[1][SA[1]], contact.impulses_dual[1][SA[1]])
	γsoc, ssoc = initialize_second_order_cone!(contact.impulses[1][SA[2;3;4]], contact.impulses_dual[1][SA[2;3;4]])
	contact.impulses[1] = [γort; γsoc]
	contact.impulses_dual[1] = [sort; ssoc]
	γort, sort = initialize_positive_orthant!(contact.impulses[2][SA[1]], contact.impulses_dual[2][SA[1]])
	γsoc, ssoc = initialize_second_order_cone!(contact.impulses[2][SA[2;3;4]], contact.impulses_dual[2][SA[2;3;4]])
	contact.impulses[2] = [γort; γsoc]
	contact.impulses_dual[2] = [sort; ssoc]
    return nothing
end

function initialize_positive_orthant!(γ::AbstractVector{T}, s::AbstractVector{T}; 
    ϵ::T = 1e-20) where T
    
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

function initialize_second_order_cone!(γ::AbstractVector{T}, s::AbstractVector{T}; 
    ϵ::T = 1e-20) where T

    e = [1.0; szeros(length(γ) - 1)] # identity element
    δs = max(-1.5 * (s[1] - norm(s[SUnitRange(2,end)])), 0)
    δγ = max(-1.5 * (γ[1] - norm(γ[SUnitRange(2,end)])), 0)

    sh = s + δs * e
    γh = γ + δγ * e
    δhs = 0.5 * transpose(sh) * γh / ((γh[1] + norm(γh[SUnitRange(2,end)])) + ϵ)
    δhγ = 0.5 * transpose(sh) * γh / ((sh[1] + norm(sh[SUnitRange(2,end)])) + ϵ)

    s0 = sh + δhs * e
    γ0 = γh + δhγ * e
	return γ0, s0
end