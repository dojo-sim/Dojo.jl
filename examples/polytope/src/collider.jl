abstract type Collider{T} end
abstract type ColliderOptions{T} end

################################################################################
# ColliderOptions
################################################################################
Dojo.@with_kw mutable struct ColliderOptions18{T} <: ColliderOptions{T}
    impact_damper::T=300.0
    impact_spring::T=100.0
    sliding_friction::T=0.2
    rolling_drag::T=0.5
    rolling_friction::T=0.04
    coulomb_smoothing::T=1000.0
    coulomb_regularizer::T=1e-3
end

################################################################################
# HalfSpaceCollider
################################################################################
mutable struct HalfSpaceCollider28{T} <: Collider{T}
    origin::Vector{T}
    normal::Vector{T}
end

function HalfSpaceCollider(origin, normal)
    normal /= norm(normal)
    return HalfSpaceCollider28(origin, normal)
end

function inside(halfspace::HalfSpaceCollider28, p)
    origin = halfspace.origin
    normal = halfspace.normal
    c = (p - origin)' * normal
    return c <= 0.0
end

################################################################################
# SoftCollider
################################################################################
mutable struct SoftCollider27{T,N} <: Collider{T}
    x::AbstractVector{T}
    q::Quaternion{T}
    mass::T
    inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    gradients::Vector{SVector{3,T}}
    collision_weights::Vector{T} # contribution to collision force F = spring_constant * collision_weight
    nerf_object::Any
    mesh::GeometryBasics.Mesh
    options::ColliderOptions{T}
end

function SoftCollider(nerf_object, mesh; N=1000, density_scale=0.1, opts=ColliderOptions18(), T=Float64)
    x = szeros(T,3)
    q = Quaternion(1,0,0,0.0)
    mass, inertia, center_of_mass = inertia_properties(nerf_object, density_scale=density_scale)
    particles, densities, gradients = sample_soft(nerf_object, N)
    collision_weigths = densities ./ sum(densities) * mass
    return SoftCollider27{T,N}(x, q, mass, inertia, center_of_mass, particles,
        densities, gradients, collision_weights, nerf_object, mesh, opts)
end

function sample_soft(nerf_object, N::Int, T=Float64, min_density=1.0, max_density=105.0)
    particles = Vector{SVector{3,T}}(undef, N)
    densities = zeros(T,N)
    gradients = Vector{SVector{3,T}}(undef, N)

    n = Int(floor((100N)^(1/3)))
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    candidate_particles = grid_particles(xrange, yrange, zrange)
    candidate_densities = py"density_query"(nerf_object, candidate_particles)

    ind = findall(x -> min_density <= x <= max_density, candidate_densities)
    good_particles = candidate_particles[ind,:]
    good_densities = candidate_densities[ind,:]
    Ng = size(good_particles)[1]
    for i = 1:N
        particles[i] = good_particles[Int(floor(i*Ng/N)),:]
        densities[i] = good_densities[Int(floor(i*Ng/N))]
        gradients[i] = finite_difference_gradient(nerf_object, particles[i])
    end
    return particles, densities, gradients
end

function finite_difference_gradient(nerf_object, particle; δ=0.01, n::Int=5)
    gradient = zeros(3)
    xrange = particle[1] .+ range(-δ, stop=δ, length=n)
    yrange = particle[2] .+ range(-δ, stop=δ, length=n)
    zrange = particle[3] .+ range(-δ, stop=δ, length=n)
    sample_particles = grid_particles(xrange, yrange, zrange)
    sample_gradients = py"density_gradient_query"(nerf_object, sample_particles)
    gradient = mean(sample_gradients, dims=1)[1,:]
    return gradient
end

################################################################################
# collision
################################################################################
function collision(halfspace::HalfSpaceCollider28{T}, collider::SoftCollider27{T,N}) where {T,N}
    ψ = 0.0
    ∇ψ = zeros(3)
    contact_normal = halfspace.normal
    barycenter = zeros(3)

    active = 0
    for i = 1:N
        particle = collider.particles[i]
        p = collider.x + Dojo.vector_rotate(particle, collider.q)
        if inside(halfspace, p)
            active +=1
            collision_weight = collider.collision_weights[i]
            gradient = collider.gradients[i]
            ψ += collision_weight
            ∇ψ += collision_weight * gradient
            barycenter += collision_weight * particle
        end
    end
    if active > 0
        ∇ψ /= ψ * N
        barycenter /= ψ
        ψ /= N
    end
    return ψ, ∇ψ, contact_normal, barycenter
end
