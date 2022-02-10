mutable struct LinearContact{T,N} <: Contact{T,N}
    cf::T
    Bx::SMatrix{4,3,T,12}
    ainv3::Adjoint{T,SVector{3,T}} # inverse matrix
    p::SVector{3,T}
    offset::SVector{3,T}

    function LinearContact(body::Body{T}, normal::AbstractVector, cf; p = szeros(T, 3), offset::AbstractVector = szeros(T, 3)) where T
        V1, V2, V3 = orthogonalcols(normal) # gives two plane vectors and the original normal axis
        A = [V1 V2 V3]
        Ainv = inv(A)
        ainv3 = Ainv[3,SA[1; 2; 3]]'
        Bx = SA{T}[
             1  0  0
            -1  0  0
             0  1  0
             0 -1  0
        ]
        new{Float64,12}(cf, Bx, ainv3, p, offset)
    end
end

function constraint(mechanism, contact::ContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:Tuple{LinearContact{T,N}}}
    bound = contact.constraints[1]
    body = get_body(mechanism, contact.parent_id)
    x2, v25, q2, ϕ25 = current_configuration_velocity(body.state)
    x3, q3 = next_configuration(body.state, mechanism.timestep)

    # transforms the velocities of the origin of the link into velocities along all 4 axes of the friction pyramid
    # vp = V(cp, B / W)_w velocity of the contact point cp, attached to body B wrt world frame, expressed in the world frame.
    vp = v25 + skew(vrotate(ϕ25, q3)) * (vrotate(bound.p, q3) - bound.offset)
    γ = contact.dual[2][1]
    sγ = contact.primal[2][1]
    ψ = contact.dual[2][2]
    sψ = contact.primal[2][2]
    β = contact.dual[2][@SVector [3,4,5,6]]
    sβ = contact.primal[2][@SVector [3,4,5,6]]
    SVector{6,T}(
        bound.ainv3 * (x3 + vrotate(bound.p,q3) - bound.offset) - sγ,
        bound.cf * γ - sum(β) - sψ,
        (bound.Bx * vp + ψ * sones(4) - sβ)...)
end

@inline function constraint_jacobian_velocity(bound::LinearContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
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

@inline function constraint_jacobian_configuration(bound::LinearContact, x3::AbstractVector, q3::UnitQuaternion,
    x2::AbstractVector, v25::AbstractVector, q2::UnitQuaternion, ϕ25::AbstractVector, λ, timestep)
    V = [bound.ainv3;
         szeros(1,3);
         szeros(4,3)]
    # Ω = FiniteDiff.finite_difference_jacobian(ϕ25 -> g(bound, s, γ, x2+timestep*v25, next_orientation(q2,ϕ25,timestep), v25, ϕ25), ϕ25)
    ∂v∂q3 = skew(vrotate(ϕ25, q3)) * ∂vrotate∂q(bound.p, q3)
    ∂v∂q3 += skew(bound.offset - vrotate(bound.p, q3)) * ∂vrotate∂q(ϕ25, q3)
    Ω = [bound.ainv3 * ∂vrotate∂q(bound.p, q3);
        szeros(1,4);
        bound.Bx * ∂v∂q3]
    return [V Ω]
end

@inline function impulse_map(bound::LinearContact, x::AbstractVector, q::UnitQuaternion, λ)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    # q * ... is a rotation by quatrnon q it is equivalent to Vmat() * Lmat(q) * Rmat(q)' * Vᵀmat() * ...
    Q = - X * q * skew(bound.p - vrotate(bound.offset, inv(q)))
    return transpose([X Q])
end

@inline function force_mapping(bound::LinearContact, x::AbstractVector, q::UnitQuaternion)
    X = [bound.ainv3;
         szeros(1,3);
         bound.Bx]
    return X
end

@inline function set_matrix_vector_entries!(mechanism::Mechanism, matrix_entry::Entry, vector_entry::Entry,
    contact::ContactConstraint{T,N,Nc,Cs,N½}) where {T,N,Nc,Cs<:Tuple{LinearContact{T,N}},N½}
    # ∇primal[dual .* primal - μ; g - s] = [diag(dual); -diag(0,1,1)]
    # ∇dual[dual .* primal - μ; g - s] = [diag(primal); -diag(1,0,0)]
    # (cf γ - ψ) dependent of ψ = dual[2][1:1]
    # B(z) * zdot - sβ dependent of sβ = primal[2][2:end]
    cf = contact.constraints[1].cf
    γ = contact.dual[2]
    s = contact.primal[2]

    ∇s1 = Diagonal(γ) # 6x6
    ∇s2 = Diagonal(-sones(T,6))
    ∇s = vcat(∇s1, ∇s2) # 12x6

    ∇γ1 = Diagonal(s) # 6x6
    ∇γ2 = @SMatrix[ 0  0  0  0  0  0;
                   cf  0 -1 -1 -1 -1;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;
                    0  1  0  0  0  0;]
    ∇γ = vcat(∇γ1, ∇γ2) # 12x6
    matrix_entry.value = hcat(∇s, ∇γ)

    # [-dual .* primal + μ; -g + s]
    vector_entry.value = vcat(-complementarityμ(mechanism, contact), -constraint(mechanism, contact))
    return
end
