mutable struct FrictionBound{T,N} <: Bound{T,N}
    cf::T

    function FrictionBound(cf, frictionid, impactid)
        new{Float64,2}(cf), frictionid, impactid
    end
end

mutable struct BetaBound{T,N} <: Bound{T,N}

    function BetaBound(frictionid)
        new{Float64,8}(), frictionid, nothing
    end
end


### Constraints and derivatives
## Position level constraints (for dynamics)
g(bound::FrictionBound, fric::Friction, Δt) = g(bound, fric.βsol[2], fric.γsolref[2])
g(bound::BetaBound, fric::Friction, Δt) = g(bound, fric.βsol[2])

# gs =  μ * γ - 1' * β - sψ
# g  =  μ * γ - 1' * β
@inline g(ffl::FrictionBound, β::AbstractVector, γ::AbstractVector) = ffl.cf*γ - SVector{1,Float64}(sum(β))

# gs = β - sη
# g  = β
@inline g(fvl::BetaBound, β::AbstractVector) = β

# contribution of the inequality constraint (impact or friction) to the friction MDP constraint d = 0 [4]
@inline function constraintForceMapping!(mechanism, fric::Friction, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{FrictionBound{T,N}}}
    # γsol is the dual of the slack variable in the friction bound
    # γsol is ψ,
    # ψ is one dimensional, that why we add .+
    # we add 1 * ψ
    # d = B(z) zdot + 1 * ψ - η
    fric.d = fric.d .+ ineqc.γsol[2]
    return
end

# contribution of the inequality constraint (impact or friction) to the friction MDP constraint d = 0 [4]
@inline function constraintForceMapping!(mechanism, fric::Friction, ineqc::InequalityConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{BetaBound{T,N}}}
    # γsol is the dual of the slack variable in the beta bound
    # γsol is η
    # we add  -η
    # d = B(z) zdot + 1 * ψ - η
    fric.d -= ineqc.γsol[2]
    return
end

## Derivatives
@inline function ∂g∂beta(::Friction{T}, ::FrictionBound) where {T}
    return [szeros(T,4) sones(T,4)]
end
@inline function ∂g∂beta(::Friction{T}, ::BetaBound) where {T}
    return [szeros(T,4,4) -I]
end

@inline function ∂gab∂ʳba(ffl::FrictionBound{T}, impact::Impact) where T
    SA{T}[0 0;0 ffl.cf], szeros(T,2,2)
end




## Complementarity
function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{FrictionBound{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2]
end

function complementarity(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{BetaBound{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2]
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{FrictionBound{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2] .- mechanism.μ
end

function complementarityμ(mechanism, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{BetaBound{T,N}},N½}
    return ineqc.γsol[2] .* ineqc.ssol[2] .- mechanism.μ
end
