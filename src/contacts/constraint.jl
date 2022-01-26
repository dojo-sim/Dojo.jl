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

function constraint_jacobian_velocity(mechanism, contact::ContactConstraint, body::Body)
    return constraint_jacobian_velocity(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

function impulse_map(mechanism, contact::ContactConstraint, body::Body)
    return impulse_map(contact.constraints[1], body, nothing, nothing, mechanism.timestep)
end

@inline function impulses!(mechanism, body::Body, contact::ContactConstraint)
    body.state.d -= impulse_map(mechanism, contact, body)' * contact.γsol[2]
    return
end

@inline function impulses_jacobian_velocity!(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    timestep = mechanism.timestep
    x3, q3 = next_configuration(body.state, timestep)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)

    for i=1:Nc
        bnd = contact.constraints[i]
        p = bnd.p
        offset = bnd.offset

        X = force_mapping(bnd)
        λ = X' * contact.γsol[2]

        ∇ = ∂pskew(VRmat(q3) * LᵀVᵀmat(q3) * λ) * -∂vrotate∂q(offset, inv(q3)) * Tmat()
        ∇ += skew(p - vrotate(offset, inv(q3))) * ∂qVRmat(LᵀVᵀmat(q3) * λ)
        ∇ += skew(p - vrotate(offset, inv(q3))) * VRmat(q3) * ∂qLᵀVᵀmat(λ)
        ∇ *= ∂integrator∂ϕ(q2, ϕ25, timestep)
        body.state.D -= [szeros(T,6,3) [szeros(T,3,3); ∇]]
    end
    return
end

@inline function off_diagonal_jacobians(mechanism, body::Body, contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs,N½}
    Z = szeros(T,N½,6)
    return [Z; -impulse_map(mechanism, contact, body)]', [Z; constraint_jacobian_velocity(mechanism, contact, body)]
end

function reset!(contact::ContactConstraint{T,N,Nc,Cs,N½}; scale::T=1.0) where {T,N,Nc,Cs,N½}
    contact.ssol[1] = scale * neutral_vector(contact.constraints[1])
    contact.ssol[2] = scale * neutral_vector(contact.constraints[1])
    contact.γsol[1] = scale * neutral_vector(contact.constraints[1])
    contact.γsol[2] = scale * neutral_vector(contact.constraints[1])
    return
end

cone_degree(contact::ContactConstraint) = sum(cone_degree.(contact.constraints))
