mutable struct InequalityConstraint{T,N,Nc,Cs,N½} <: Constraint{T,N}
    id::Int64
    name::String

    # Currently only single constraint and child
    constraints::Cs # can be of type
    parentid::Int64
    childids::SVector{1,Union{Int64,Nothing}}
    ssol::Vector{SVector{N½,T}} # holds the slack variable
    γsol::Vector{SVector{N½,T}} # holds the dual of the slack variable

    function InequalityConstraint(data; name::String="")
        bound, parentid, childid = data
        T = getT(bound)

        childids = [childid]
        constraint = Tuple([bound])
        N = length(constraint[1])
        N½ = Int64(N/2)

        ssol = [neutral_vector(bound) for i = 1:2]
        γsol = [neutral_vector(bound) for i = 1:2]
        new{T,N,1,typeof(constraint),N½}(getGlobalID(), name, constraint, parentid, childids, ssol, γsol)
    end
end

function resetVars!(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    ineqc.ssol[1] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.ssol[2] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.γsol[1] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.γsol[2] = scale * neutral_vector(ineqc.constraints[1])
    return
end

# contribution of the inequality constraint (impact or friction) to the dynamics equation d = 0
@inline function impulses!(mechanism, body::Body, ineqc::InequalityConstraint)
    body.state.d -= ∂g∂ʳpos(mechanism, ineqc, body)' * ineqc.γsol[2]
    return
end

@inline function ∂impulses∂v!(mechanism, body::Body, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Δt = mechanism.Δt
    x3, q3 = posargs3(body.state, Δt)
    x2, v25, q2, ϕ25 = fullargssol(body.state)

    for i=1:Nc
        bnd = ineqc.constraints[i]
        p = bnd.p
        offset = bnd.offset

        X = forcemapping(bnd)
        λ = X' * ineqc.γsol[2]

        ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
        ∇ += skew(p - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
        ∇ += skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
        ∇ *= ∂integrator∂ϕ(q2, ϕ25, Δt)
        body.state.D -= [szeros(T,6,3) [szeros(T,3,3); ∇]]
    end
    return
end

@inline function ∂gab∂ʳba(mechanism, body::Body, ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; -∂g∂ʳpos(mechanism, ineqc, body)]', [Z; ∂g∂ʳvel(mechanism, ineqc, body)]
end
@inline function ∂gab∂ʳba(mechanism, ineqc1::InequalityConstraint, ineqc2::InequalityConstraint)
    G1, G2 = ∂gab∂ʳba(ineqc1.constraints[1], ineqc2.constraints[1])
    return G1, G2
end

function ∂g∂ʳposa(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳposa(ineqc.constraints[1], body, nothing, nothing, mechanism.Δt)
end

function ∂g∂ʳvela(mechanism, ineqc::InequalityConstraint, body::Body)
    return ∂g∂ʳvela(ineqc.constraints[1], body, nothing, nothing, mechanism.Δt)
end

function cone_degree(ineqc::InequalityConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    return sum(cone_degree.(ineqc.constraints))
end
