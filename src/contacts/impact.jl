mutable struct ImpactContact{T,N} <: Contact{T,N}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function ImpactContact(body::Body{T}, normal::AbstractVector; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(ainv3, p, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{ImpactContact{T,N}}}
    bound = contact.constraints[1]
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    SVector{1,T}(bound.ainv3 * (x3 + vrotate(bound.p,q3) - bound.offset) - contact.primal[2][1])
end

@inline function constraint_jacobian_velocity(bound::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    V = bound.ainv3 * timestep
    Ω = bound.ainv3 * ∂vrotate∂q(bound.p, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    return [V Ω]
end

@inline function constraint_jacobian_configuration(bound::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    X = bound.ainv3
    Q = bound.ainv3 * ∂vrotate∂q(bound.p, q3)
    return [X Q]
end

@inline function impulse_map(bound::ImpactContact, x::AbstractVector, q::UnitQuaternion, λ)
    X = bound.ainv3
    # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return transpose([X Q])
end

@inline function force_mapping(bound::ImpactContact, x::AbstractVector, q::UnitQuaternion)
    X = bound.ainv3
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{ImpactContact{T,N}},N½}
    # ∇primal[dual .* primal - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* primal - μ; g - s] = [diag(primal); -diag(1,0,0)]
    γ = contact.dual[2]
    s = contact.primal[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* primal + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end
