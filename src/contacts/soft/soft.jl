################################################################################
# ColliderOptions
################################################################################
abstract type ColliderOptions{T} end
@with_kw mutable struct ColliderOptions250{T} <: ColliderOptions{T}
    impact_damper::T=1e7
    impact_spring::T=1e7
    sliding_drag::T=0.0
    sliding_friction::T=0.1
    rolling_drag::T=0.05
    rolling_friction::T=0.01
    coulomb_smoothing::T=1000.0
    coulomb_regularizer::T=1e-3
end

abstract type Collider{T} end
mutable struct SoftCollider250{T,N} <: Collider{T}
    # x::AbstractVector{T}
    # q::Quaternion{T}
    # mass::T
    # inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    density_gradients::Vector{SVector{3,T}}
    weights::Vector{T} # contribution to collision force F = spring_constant * collision_weight
    weight_gradients::Vector{SVector{3,T}}
    nerf_object::Any
    # mesh::GeometryBasics.Mesh
    options::ColliderOptions{T}
end

function SoftCollider(nerf_object; N=1000, density_scale=0.1, opts=ColliderOptions250(), T=Float64)
    # x = szeros(T,3)
    # q = Quaternion(1,0,0,0.0)
    mass, inertia, center_of_mass = inertia_properties(nerf_object, density_scale=density_scale)
    particles, densities, density_gradients = sample_soft(nerf_object, N)
    weights = densities ./ sum(densities) * mass
    weight_gradients = density_gradients ./ sum(densities) * mass
    return SoftCollider250{T,N}(
        # x, q, mass, inertia,
        center_of_mass,
        particles, densities, density_gradients,
        weights, weight_gradients, nerf_object,
        # mesh,
        opts)
end

################################################################################
# HalfSpaceCollider
################################################################################
mutable struct HalfSpaceCollider250{T} <: Collider{T}
    origin::SVector{3,T}
    normal::SVector{3,T}
end

function HalfSpaceCollider(origin, normal)
    normal /= norm(normal)
    return HalfSpaceCollider250(origin, normal)
end

function inside(collider::HalfSpaceCollider250, p)
    origin = collider.origin
    normal = collider.normal
    c = (p - origin)' * normal
    return c <= 0.0
end

function halfspace_collision(collider::SoftCollider250{T,N}, xp, qp, xc, qc) where {T,N}
    center_of_mass = collider.center_of_mass
    Ψ = Vector{T}()
    active_particles = []
    barycenter = szeros(T,3)

    halfspace_origin = xc
    halfspace_normal = vector_rotate([0,0,1.0], qc)
    halfspace = HalfSpaceCollider(halfspace_origin, halfspace_normal)

    for i = 1:N
        particle = collider.particles[i]
        p = xp + Dojo.vector_rotate(particle - center_of_mass, qp)
        if inside(halfspace, p)
            push!(Ψ, collider.weights[i])
            push!(active_particles, particle)
        end
    end

    num_active = length(Ψ)
    ψ = sum(Ψ)
    (ψ > 0) && (barycenter = sum(Vector{SVector{3,T}}([Ψ[i] * active_particles[i] for i=1:num_active])) / ψ)
    contact_normal = collision_normal(halfspace, collider, barycenter)
    return ψ, barycenter, contact_normal
end

function collision_normal(collider::HalfSpaceCollider250{T}, soft_collider::SoftCollider250{T,N}, particle) where {T,N}
    collider.normal
end











"""
    SoftContact{T,N} <: Contact{T,N}

    contact object for soft contact

    collision: Collision
"""
abstract type SoftContact{T,N} <: Contact{T,N} end
mutable struct Soft250Contact{T,N} <: SoftContact{T,N}
    friction_parameterization::SMatrix{0,2,T,0}
    collision::Collision{T,0,3,0}
    contact_origin::SVector{3,T} # origin of the frame in which the nerf is defined
    collider::Collider{T}
end

function SoftContact(body::Body{T}, normal::AbstractVector{T};
    collider::Collider=SOFT,
    contact_origin=szeros(T, 3),
    contact_radius=0.0) where T

    # contact directions
    V1, V2, V3 = orthogonal_columns(normal) #
    A = [V1 V2 V3]
    Ainv = inv(A)
    contact_normal = Ainv[3, SA[1; 2; 3]]'

    # friction parametrization
    parameterization = szeros(T, 0, 2)

    # collision
    collision = SphereHalfSpaceCollision(szeros(T, 0, 3), contact_normal, SVector{3}(contact_origin), contact_radius)

    Soft250Contact{Float64,1}(parameterization, collision, contact_origin, collider)
end

function constraint(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    # contact model
    model = contact.model
    collider = model.collider

    # parent
    pbody = get_body(mechanism, contact.parent_id)
    x2p, v25p, q2p, ϕ25p = current_configuration_velocity(pbody.state)

    # child
    cbody = get_body(mechanism, contact.child_id)
    x2c, q2c = current_configuration(cbody.state)

    # collision
    ψ, barycenter, contact_normal = halfspace_collision(collider, x2p, q2p, x2c, q2c)
    # velocities
    vc = v25p + vector_rotate(skew(collider.center_of_mass - barycenter) * ϕ25p, q2p) # TODO could be q3
    vc_normal = vc' * contact_normal# * contact_normal
    # vc_tangential = vc - vc_normal
    @show ψ
    # constraint
    # SVector{1,T}(distance(model.collision, x3p, q3p, x3c, q3c) - contact.impulses_dual[2][1])
    SVector{1,T}(ψ * collider.options.impact_spring * vc_normal - contact.impulses[2][1])
end

function constraint_jacobian(contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    ∇γ = -Diagonal(sones(N))
    return ∇γ
end

function constraint_jacobian_configuration(relative::Symbol, model::SoftContact,
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where T
    return szeros(T,1,6)
end

function constraint_jacobian_velocity(relative::Symbol, model::SoftContact{T},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where T

    collider = model.collider
    impact_spring = collider.options.impact_spring
    x2p = next_position(xp, -vp, timestep)
    q2p = next_orientation(qp, -ϕp, timestep)
    x2c = next_position(xc, -vc, timestep)
    q2c = next_orientation(qc, -ϕc, timestep)

    ψ, barycenter, contact_normal = halfspace_collision(collider, x2p, q2p, x2c, q2c)

    # recover current orientation
    if relative == :parent
        V = ψ * impact_spring * contact_normal'
        Ω = ψ * impact_spring * contact_normal' * rotation_matrix(q2p) * skew(collider.center_of_mass - barycenter)
    elseif relative == :child
        V = szeros(T,1,3)
        Ω = szeros(T,1,3)
    end
    return [V Ω]
end

function force_mapping(relative::Symbol, model::SoftContact,
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion)

    # X = contact_normal(model.collision, xp, qp, xc, qc)'
    X = [0, 0, 1.0] # TODO hard-coded
    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end

function ∂force_mapping_jvp∂x(relative::Symbol, jacobian::Symbol,
    model::SoftContact{T},
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector) where T

    # X = ∂contact_normal_transpose∂x(jacobian, model.collision, xp, qp, xc, qc) * λ[1]
    X = szeros(T,3,3)
    @show size(X)
    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end

function ∂force_mapping_jvp∂q(relative::Symbol, jacobian::Symbol,
    model::SoftContact{T},
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion,
    λ::AbstractVector) where T

    # X = ∂contact_normal_transpose∂q(jacobian, model.collision, xp, qp, xc, qc) * λ[1]
    X = szeros(T,3,4)

    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end
