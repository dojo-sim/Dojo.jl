mutable struct ImpactContact{T,N} <: Contact{T,N}
    surface_normal_projector::Adjoint{T,SVector{3,T}} # inverse matrix
    contact_point::SVector{3,T}
    offset::SVector{3,T}

    function ImpactContact(body::Body{T}, normal::AbstractVector; contact_point = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        surface_normal_projector = Ainv[3,SA[1; 2; 3]]'
        new{Float64,2}(surface_normal_projector, contact_point, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:ImpactContact{T,N}}
    model = contact.model
    body = get_body(mechanism, contact.parent_id)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    SVector{1,T}(model.surface_normal_projector * (x3 + vrotate(model.contact_point,q3) - model.offset) - contact.impulses_dual[2][1])
end

@inline function constraint_jacobian_velocity(model::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    V = model.surface_normal_projector * timestep
    Ω = model.surface_normal_projector * ∂vrotate∂q(model.contact_point, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep)
    return [V Ω]
end

@inline function constraint_jacobian_configuration(model::ImpactContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    X = model.surface_normal_projector
    Q = model.surface_normal_projector * ∂vrotate∂q(model.contact_point, q3)
    return [X Q]
end

@inline function force_mapping(model::ImpactContact, x::AbstractVector, q::UnitQuaternion)
    X = model.surface_normal_projector
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:ImpactContact{T,N},N½}
    # ∇impulses[dual .* impulses - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* impulses - μ; g - s] = [diag(impulses); -diag(1,0,0)]
    γ = contact.impulses[2]
    s = contact.impulses_dual[2]

    ∇s = hcat(γ, -Diagonal(sones(N½)))
    ∇γ = hcat(s, -Diagonal(szeros(N½)))
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* impulses + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end
