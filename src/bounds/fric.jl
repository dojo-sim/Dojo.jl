



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











mutable struct Friction{T} <: Component{T}
    id::Int64
    name::String

    Bx::SMatrix{4,3,T,12}
    p::SVector{3,T}

    parentid::Int64
    childids::SVector{2,Int64}

    βsol::Vector{SVector{4,T}}
    γsolref::Vector{SVector{1,T}} #TODO this is a reference to the associated impact γ to avoid allocations

    d::SVector{4,T}

    function Friction(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3), name::String="") where T
        N = 4

        Bx = SA{T}[
            1 0 0
            -1 0 0
            0 1 0
            0 -1 0
        ]

        frictionid = getGlobalID()
        impact = InequalityConstraint(Impact(body, normal; p = p, offset = offset))
        frictionbound = InequalityConstraint(FrictionBound(cf, frictionid, impact.id))
        betabound = InequalityConstraint(BetaBound(frictionid))

        βsol = [szeros(T, N) for i=1:2]
        γsolref = impact.γsol

        d = szeros(T, 4)

        new{T}(frictionid, name, Bx, p, body.id, [frictionbound.id; betabound.id], βsol, γsolref, d), [impact; frictionbound; betabound]
    end
end


Base.length(::Friction) = 4

@inline ∂g∂ʳself(mechanism, fric::Friction{T}) where {T} = szeros(T,4,4)

@inline Bq(Bxmat, p, q) = Bxmat*VRᵀmat(q)*LVᵀmat(q)*skew(-p)
@inline Bmat(Bxmat, p, q) = [Bxmat Bq(Bxmat, p, q)]

function g(mechanism, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    x, v, q, ω = fullargssol(body.state)
    Bxmat = fric.Bx
    Bqmat = Bq(Bxmat, fric.p, q)

    fric.d = Bxmat*v + Bqmat*ω

    for childid in fric.childids
        constraintForceMapping!(mechanism, fric, getineqconstraint(mechanism, childid))
    end

    return fric.d
end

@inline function constraintForceMapping!(mechanism, body::Body, fric::Friction)
    body = getbody(mechanism, fric.parentid)
    x, q = posargsk(body.state)

    body.state.d -= Bmat(fric.Bx, fric.p, q)' * fric.βsol[2]
    return
end

@inline function ∂gab∂ʳba(mechanism, body::Body, fric::Friction)
    x, q = posargsk(body.state)
    B = Bmat(fric.Bx, fric.p, q)

    return -B', B
end

@inline function ∂gab∂ʳba(mechanism, fric::Friction, ineqc::InequalityConstraint)
    G = ∂g∂beta(fric, ineqc.constraints[1])

    return G, -G'
end











mutable struct Impact{T,N} <: Bound{T,N}
    ainv3::Adjoint{T,SVector{3,T}}
    p::SVector{3,T}
    offset::SVector{3,T}


    function Impact(body::Body{T}, normal::AbstractVector; p::AbstractVector = zeros(3),  offset::AbstractVector = zeros(3)) where T
        # Derived from plane equation a*v1 + b*v2 + distance*v3 = p - offset
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'

        new{T,2}(ainv3, p, offset), body.id, nothing
    end
end


### Constraints and derivatives
## Position level constraints (for dynamics)
@inline g(impact::Impact{T}, x::AbstractVector, q::UnitQuaternion) where T = SVector{1,T}(impact.ainv3 * (x + vrotate(impact.p,q) - impact.offset))

## Derivatives NOT accounting for quaternion specialness
@inline function ∂g∂pos(impact::Impact, x::AbstractVector, q::UnitQuaternion)
    p = impact.p
    X = impact.ainv3
    Q = impact.ainv3 * (VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p)))
    return X, Q
end
