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

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:ImpactContact{T,N}}
    model = contact.model
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    SVector{1,T}(model.ainv3 * (x3 + vrotate(model.p,q3) - model.offset) - contact.impulses[2][1])
end

@inline function constraint_jacobian_velocity(model::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    V = model.ainv3 * timestep
    Ω = model.ainv3 * ∂vrotate∂q(model.p, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    return [V Ω]
end

@inline function constraint_jacobian_configuration(model::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    X = model.ainv3
    Q = model.ainv3 * ∂vrotate∂q(model.p, q3)
    return [X Q]
end

@inline function impulse_map(model::ImpactContact, x::AbstractVector, q::UnitQuaternion, λ)
    X = model.ainv3
    # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(model.p - vrotate(model.offset, inv(q)))
    return [X Q]
end

@inline function force_mapping(model::ImpactContact)
    X = model.ainv3
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:ImpactContact{T,N},N½}
    # ∇impulses[dual .* impulses - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* impulses - μ; g - s] = [diag(impulses); -diag(1,0,0)]
    γ = contact.impulses_dual[2]
    s = contact.impulses[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* impulses + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end
