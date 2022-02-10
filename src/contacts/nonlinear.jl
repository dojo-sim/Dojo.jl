mutable struct NonlinearContact{T,N} <: Contact{T,N}
    cf::T
    Bx::SMatrix{2,3,T,6}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function NonlinearContact(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        Bx = SA{T}[
            1 0 0
            0 1 0
        ]
        new{T,8}(cf, Bx, ainv3, p, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}}}
    bound = contact.constraints[1]
    body = get_body(mechanism, contact.parent_id)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, mechanism.timestep)
    constraint(bound, contact.primal[2], contact.dual[2], x3, q3, v25, ϕ25)
end

function constraint(bound::NonlinearContact, s::AbstractVector{T}, γ::AbstractVector{T},
        x3::AbstractVector{T}, q3::UnitQuaternion{T}, v25::AbstractVector{T},
        ϕ25::AbstractVector{T}) where T

    # transforms the velocities of the origin of the link into velocities
    vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
    SVector{4,T}(
        bound.ainv3 * (x3 + vrotate(bound.p, q3) - bound.offset) - s[1],
        bound.cf * γ[1] - γ[2],
        (bound.Bx * vp - s[@SVector [3,4]])...)
end

@inline function constraint_jacobian_velocity(bound::NonlinearContact{T}, x3::AbstractVector{T}, q3::UnitQuaternion{T},
    x2::AbstractVector{T}, v25::AbstractVector{T}, q2::UnitQuaternion{T}, ϕ25::AbstractVector{T}, λ, timestep::T) where T
    V = [bound.ainv3 * timestep;
        szeros(1,3);
        bound.Bx]
    # Ω = FiniteDiff.finite_difference_jacobian(ϕ25 -> g(bound, s, γ, x2+timestep*v25, next_orientation(q2,ϕ25,timestep), v25, ϕ25), ϕ25)
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    ∂v∂ϕ25 = skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂p(ϕ25, q3)
    Ω = [bound.ainv3 * ∂vrotate∂q(bound.p, q3) * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep);
        szeros(1,3);
        bound.Bx * (∂v∂ϕ25 + ∂v∂q3 * rotational_integrator_jacobian_velocity(q2, ϕ25, timestep))]
    return [V Ω]
end

@inline function constraint_jacobian_configuration(bound::NonlinearContact{T}, x3::AbstractVector{T}, q3::UnitQuaternion{T},
    x2::AbstractVector{T}, v25::AbstractVector{T}, q2::UnitQuaternion{T}, ϕ25::AbstractVector{T}, λ, timestep::T) where T
    X = [bound.ainv3;
        szeros(1,3);
        szeros(2,3)]
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    Q = [bound.ainv3 * ∂vrotate∂q(bound.p, q3);
        szeros(1,4);
        bound.Bx * ∂v∂q3]
    return [X Q]
end

@inline function impulse_map(bound::NonlinearContact, x::AbstractVector, q::UnitQuaternion, λ)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    # q * ... is a rotation by quaternion q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return transpose([X Q])
end

@inline function force_mapping(bound::NonlinearContact, x::AbstractVector, q::UnitQuaternion)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    # ∇primal[dual .* primal - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* primal - μ; g - s] = [diag(primal); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = dual[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = primal[2][2:end]
    cf = contact.constraints[1].cf
    γ = contact.dual[2] + 1e-10*neutral_vector(contact.constraints[1]) # TODO need to check this is legit
    s = contact.primal[2] + 1e-10*neutral_vector(contact.constraints[1]) # TODO need to check this is legit

    # ∇s = [contact.dual[2][1] szeros(1,3); szeros(3,1) cone_product_jacobian(contact.dual[2][2:4]); Diagonal([-1, 0, -1, -1])]
    ∇s1 = [γ[SA[1]]; szeros(T,3)]'
    ∇s2 = [szeros(T,3,1) cone_product_jacobian(γ[SA[2,3,4]])]
    ∇s3 = Diagonal(SVector{4,T}(-1, 0, -1, -1))
    ∇s = [∇s1; ∇s2; ∇s3]

    # ∇γ = [contact.primal[2][1] szeros(1,3); szeros(3,1) cone_product_jacobian(contact.primal[2][2:4]); szeros(1,4); cf -1 0 0; szeros(2,4)]
    ∇γ1 = [s[SA[1]]; szeros(T,3)]'
    ∇γ2 = [szeros(T,3,1) cone_product_jacobian(s[SA[2,3,4]])]
    ∇γ3 = SA[0   0 0 0;
             cf -1 0 0;
             0   0 0 0;
             0   0 0 0;]
    ∇γ = [∇γ1; ∇γ2; ∇γ3]

    # matrix_entry.value = [[contact.dual[2][1] szeros(1,3); szeros(3,1) cone_product_jacobian(contact.dual[2][2:4]); Diagonal([-1, 0, -1, -1])] [contact.primal[2][1] szeros(1,3); szeros(3,1) cone_product_jacobian(contact.primal[2][2:4]); szeros(1,4); cf -1 0 0; szeros(2,4)]]
    matrix_entry.value = [∇s ∇γ]

    # [-dual .* primal + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end

function complementarity(mechanism, contact::ContactConstraint{T,N,Nc,Cs,N½};
        scaling::Bool = false) where {T,N,Nc,Cs<:Tuple{NonlinearContact{T,N}},N½}
    γ = contact.dual[2]
    s = contact.primal[2]
    return vcat(γ[1] * s[1], cone_product(γ[@SVector [2,3,4]], s[@SVector [2,3,4]]))
end

neutral_vector(bound::NonlinearContact{T,N}) where {T,N} = [sones(T, 2); szeros(T, Int(N/2) -2)]

cone_degree(bound::NonlinearContact) = 2
