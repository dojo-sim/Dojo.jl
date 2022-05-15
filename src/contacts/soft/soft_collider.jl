################################################################################
# ColliderOptions
################################################################################
@with_kw mutable struct ColliderOptions{T}
    impact_damper::T=3e5
    impact_spring::T=3e4
    sliding_drag::T=0.1
    sliding_friction::T=0.2
    rolling_drag::T=0.0
    rolling_friction::T=0.01
    coulomb_smoothing::T=3e1
    coulomb_regularizer::T=1e-3
end

abstract type Collider{T} end
mutable struct SoftCollider{T,N} <: Collider{T}
    mass::T
    inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    density_gradients::Vector{SVector{3,T}}
    weights::Vector{T} # contribution to collision force F = spring_constant * collision_weight
    weight_gradients::Vector{SVector{3,T}}
    nerf_object::Any
    normalizer
    options::ColliderOptions{T}
end

function SoftCollider(nerf_object, normalizer; num_particle=1000,
        opts=ColliderOptions(), T=Float64)

    mass, inertia, center_of_mass = inertia_properties(nerf_object,
        normalizer=normalizer)
    particles, densities, density_gradients = sample_soft(nerf_object, num_particle,
        normalizer=normalizer, particle_noise=0.005)

    weights = densities ./ sum(densities) * mass
    weight_gradients = density_gradients ./ sum(densities) * mass

    return SoftCollider{T,num_particle}(
        mass,
        inertia,
        center_of_mass,
        particles,
        densities,
        density_gradients,
        weights,
        weight_gradients,
        nerf_object,
        normalizer,
        opts,
        )
end


function SoftCollider(nerf::Symbol; opts=ColliderOptions(), T=Float64)

    dir = joinpath(osf_loader_dir(), "assets/collider", String(nerf) * ".jld2")
    file = JLD2.jldopen(dir)
    mass = file["mass"]
    inertia = file["inertia"]
    center_of_mass = file["center_of_mass"]
    particles = file["particles"]
    densities = file["densities"]
    density_gradients = file["density_gradients"]
    JLD2.close(file)

    num_particle = size(particles)[1]
    weights = densities ./ sum(densities) * mass
    weight_gradients = density_gradients ./ sum(densities) * mass

    return SoftCollider{T,num_particle}(
        mass,
        inertia,
        center_of_mass,
        particles,
        densities,
        density_gradients,
        weights,
        weight_gradients,
        nerf_object,
        normalizer,
        opts,
        )
end
