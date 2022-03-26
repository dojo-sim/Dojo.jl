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
    coulomb_smoothing::T=10.0
    coulomb_regularizer::T=1e-3
end

abstract type Collider{T} end
mutable struct SoftCollider250{T,N} <: Collider{T}
    # mass::T
    # inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    density_gradients::Vector{SVector{3,T}}
    weights::Vector{T} # contribution to collision force F = spring_constant * collision_weight
    weight_gradients::Vector{SVector{3,T}}
    nerf_object::Any
    options::ColliderOptions{T}
end

function SoftCollider(nerf_object; N=1000, density_scale=0.1, opts=ColliderOptions250(), T=Float64)
    mass, inertia, center_of_mass = inertia_properties(nerf_object, density_scale=density_scale)
    particles, densities, density_gradients = sample_soft(nerf_object, N, particle_noise=0.005)
    weights = densities ./ sum(densities) * mass
    weight_gradients = density_gradients ./ sum(densities) * mass
    return SoftCollider250{T,N}(
        # mass, inertia,
        center_of_mass,
        particles, densities, density_gradients,
        weights, weight_gradients, nerf_object,
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
    VERBOSE && num_active
    ψ = sum(Ψ)
    (ψ > 0) && (barycenter = sum(Vector{SVector{3,T}}([Ψ[i] * active_particles[i] for i=1:num_active])) / ψ)
    contact_normal = collision_normal(halfspace, collider, barycenter)
    return ψ, barycenter, contact_normal
end

function collision_normal(collider::HalfSpaceCollider250{T}, soft_collider::SoftCollider250{T,N}, particle) where {T,N}
    collider.normal
end

coulomb_direction(v, smoothing=1e3, regularizer=1e-3) = - atan(smoothing * norm(v)) * v/(regularizer + norm(v))
function ∂coulomb_direction∂v(v, smoothing=1e3, regularizer=1e-3)
    ∇ = - 1 / (1 + smoothing^2 * v'*v) * smoothing * v/(regularizer + norm(v)) * v'/(norm(v)+1e-20) +
        - atan(smoothing * norm(v)) * (1/(norm(v) + regularizer) * Diagonal(sones(3)) - v*v' ./ ((norm(v)+1e-20) * (norm(v) + regularizer)^2))
    return ∇
end

v = srand(3)
coulomb_direction(v)
J0 = ∂coulomb_direction∂v(v)
J1 = FiniteDiff.finite_difference_jacobian(v-> coulomb_direction(v), v)
norm(J0 - J1, Inf) < 1e-5









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
    contact_origin=-SOFT.center_of_mass,
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

    Soft250Contact{Float64,6}(parameterization, collision, contact_origin, collider)
end

function constraint(mechanism, contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    timestep = mechanism.timestep
    pbody = get_body(mechanism, contact.parent_id)
    cbody = get_body(mechanism, contact.child_id)
    xp, vp, qp, ϕp = next_configuration_velocity(pbody.state, timestep)
    xc, vc, qc, ϕc = next_configuration_velocity(cbody.state, timestep)
    impulse = constraint(contact.model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)

    VERBOSE && println(
        # "x2p:", scn.(x2p),
        # "    ψ:", scn(ψ),
        "    λ:", scn.(contact.impulses[2]),
        "    imp:", scn.(impulse),
        # "    err:", scn(impulse[1] - contact.impulses[2][1])
        )
    # @show scn.(impulse)
    # @show scn.(contact.impulses[2])
    # @show scn.(impulse - contact.impulses[2])
    return impulse - contact.impulses[2]
end

function constraint(model::SoftContact{T,N},
        xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
        xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
        timestep; intermediate::Bool=false) where {T,N}

    collider = model.collider
    opts = collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    x2p = next_position(xp, -vp, timestep)
    q2p = next_orientation(qp, -ϕp, timestep)
    x2c = next_position(xc, -vc, timestep)
    q2c = next_orientation(qc, -ϕc, timestep)

    # collision
    ψ, barycenter, contact_normal = halfspace_collision(collider, x2p, q2p, x2c, q2c)
    # velocities
    vc = vp + vector_rotate(skew(collider.center_of_mass - barycenter) * ϕp, q2p) # TODO could be q3p
    vc_normal = contact_normal * contact_normal' * vc
    vc_tangential = vc - vc_normal
    # impact
    F_impact = ψ * opts.impact_spring * contact_normal
    F_impact -= ψ * opts.impact_damper * vc_normal
    # friction
    F_friction = -opts.sliding_drag * norm(F_impact) * vc_tangential
    F_friction += opts.sliding_friction * norm(F_impact) * coulomb_direction(vc_tangential, smoothing, regularizer)
    F_contact = F_impact + F_friction

    # τ2 = Dojo.skew(barycenter - center_of_mass) * Dojo.vector_rotate(F_contact, inv(q2)) # TODO could be q3p
    τ_contact = -opts.rolling_drag * norm(F_impact) * ϕp
    τ_contact += opts.rolling_friction * norm(F_impact) * coulomb_direction(ϕp, smoothing, regularizer)

    # constraint
    impulse = timestep * [F_contact; τ_contact]
    !intermediate && return impulse
    return ψ, barycenter, contact_normal, vc_normal, vc_tangential, F_impact, F_friction, impulse
end

function constraint_jacobian(contact::SoftContactConstraint{T,N,Nc,Cs}) where {T,N,Nc,Cs<:SoftContact{T,N}}
    ∇γ = -Diagonal(sones(N))
    return ∇γ
end

function constraint_jacobian_configuration(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}
    return szeros(T,N,6)
end

function constraint_jacobian_velocity(relative::Symbol, model::SoftContact{T,N},
    xp::AbstractVector, vp::AbstractVector, qp::Quaternion, ϕp::AbstractVector,
    xc::AbstractVector, vc::AbstractVector, qc::Quaternion, ϕc::AbstractVector,
    timestep) where {T,N}

    collider = model.collider
    opts = model.collider.options
    smoothing = opts.coulomb_smoothing
    regularizer = opts.coulomb_regularizer
    ψ, barycenter, contact_normal, vc_normal, vc_tangential, F_impact, F_friction, impulse = constraint(model,
        xp, vp, qp, ϕp,
        xc, vc, qc, ϕc,
        timestep, intermediate=true)
    # recover current orientation
    q2p = next_orientation(qp, -ϕp, timestep)
    ∂vc∂v = Diagonal(sones(3))
    ∂vc∂ϕ = rotation_matrix(q2p) * skew(collider.center_of_mass - barycenter)
    if relative == :parent
        # impact
        ∂F_impact∂vc = -ψ * opts.impact_damper * contact_normal * contact_normal'
        FV_impact = ∂F_impact∂vc * ∂vc∂v
        FΩ_impact = ∂F_impact∂vc * ∂vc∂ϕ
        # friction
        ∂F_friction∂vc =  -opts.sliding_drag * vc_tangential * (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc += -opts.sliding_drag * norm(F_impact) * (Diagonal(sones(3)) - contact_normal * contact_normal')
        ∂F_friction∂vc +=  opts.sliding_friction *
            coulomb_direction(vc_tangential, smoothing, regularizer) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        ∂F_friction∂vc +=  opts.sliding_friction * norm(F_impact) *
            ∂coulomb_direction∂v(vc_tangential, smoothing, regularizer) *
            (Diagonal(sones(3)) - contact_normal * contact_normal')
        FV_friction = ∂F_friction∂vc * ∂vc∂v
        FΩ_friction = ∂F_friction∂vc * ∂vc∂ϕ

        τV_contact = (-opts.rolling_drag * ϕp + opts.rolling_friction * coulomb_direction(ϕp, smoothing, regularizer)) *
            (FV_impact * F_impact/(norm(F_impact) + 1e-20))'
        τΩ_contact = τV_contact * ∂vc∂ϕ
        τΩ_contact += norm(F_impact) * (-opts.rolling_drag * Diagonal(sones(3)) +
            opts.rolling_friction * ∂coulomb_direction∂v(ϕp, smoothing, regularizer))
        # @show ∂coulomb_direction∂v(ϕp, smoothing, regularizer)
        # contact
        V = [FV_impact + FV_friction; τV_contact]
        Ω = [FΩ_impact + FΩ_friction; τΩ_contact]
        # @show scn.(V)
        # @show scn.(Ω)
    elseif relative == :child
        V = szeros(T,N,3)
        Ω = szeros(T,N,3)
    end
    return timestep * [V Ω]
end

function force_mapping(relative::Symbol, model::SoftContact,
    xp::AbstractVector, qp::Quaternion,
    xc::AbstractVector, qc::Quaternion)

    X = Diagonal(sones(3)) # TODO hard-coded
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

    X = szeros(T,3,3)
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

    X = szeros(T,3,4)

    if relative == :parent
        return X
    elseif relative == :child
        return -1.0 * X
    end
end
