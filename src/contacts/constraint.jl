mutable struct ContactConstraint{T,N,Nc,Cs,N½} <: Constraint{T,N}
    id::Int64
    name::Symbol

    # Currently only single constraint and child
    constraints::Cs # can be of type
    parentid::Int64
    childids::SVector{1,Union{Int64,Nothing}}
    ssol::Vector{SVector{N½,T}} # holds the slack variable
    γsol::Vector{SVector{N½,T}} # holds the dual of the slack variable

    function ContactConstraint(data; name::Symbol=Symbol("contact_" * randstring(4)))
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

function ∂g∂v(mechanism, ineqc::ContactConstraint, body::Body)
    return ∂g∂v(ineqc.constraints[1], body, nothing, nothing, mechanism.Δt)
end

function G(mechanism, ineqc::ContactConstraint, body::Body)
    return G(ineqc.constraints[1], body, nothing, nothing, mechanism.Δt)
end

@inline function impulses!(mechanism, body::Body, ineqc::ContactConstraint)
    body.state.d -= G(mechanism, ineqc, body)' * ineqc.γsol[2]
    return
end

@inline function ∂impulses∂v!(mechanism, body::Body, ineqc::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Δt = mechanism.Δt
    x3, q3 = next_configuration(body.state, Δt)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)

    for i=1:Nc
        bnd = ineqc.constraints[i]
        p = bnd.p
        offset = bnd.offset

        X = force_mapping(bnd)
        λ = X' * ineqc.γsol[2]

        ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
        ∇ += skew(p - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
        ∇ += skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
        ∇ *= ∂integrator∂ϕ(q2, ϕ25, Δt)
        body.state.D -= [szeros(T,6,3) [szeros(T,3,3); ∇]]
    end
    return
end

@inline function ∂gab∂ʳba(mechanism, body::Body, ineqc::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; -G(mechanism, ineqc, body)]', [Z; ∂g∂v(mechanism, ineqc, body)]
end

function reset!(ineqc::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    ineqc.ssol[1] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.ssol[2] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.γsol[1] = scale * neutral_vector(ineqc.constraints[1])
    ineqc.γsol[2] = scale * neutral_vector(ineqc.constraints[1])
    return
end

cone_degree(ineqc::ContactConstraint) = sum(cone_degree.(ineqc.constraints))
