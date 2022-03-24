abstract type Collider{T} end
abstract type ColliderOptions{T} end

################################################################################
# ColliderOptions
################################################################################
@with_kw mutable struct ColliderOptions110{T} <: ColliderOptions{T}
    impact_damper::T=1e7
    impact_spring::T=1e7
    sliding_drag::T=0.0
    sliding_friction::T=0.1
    rolling_drag::T=0.05
    rolling_friction::T=0.01
    coulomb_smoothing::T=1000.0
    coulomb_regularizer::T=1e-3
end

################################################################################
# HalfSpaceCollider
################################################################################
mutable struct HalfSpaceCollider110{T} <: Collider{T}
    origin::Vector{T}
    normal::Vector{T}
end

function HalfSpaceCollider(origin, normal)
    normal /= norm(normal)
    return HalfSpaceCollider110(origin, normal)
end

function inside(collider::HalfSpaceCollider110, p)
    origin = collider.origin
    normal = collider.normal
    c = (p - origin)' * normal
    return c <= 0.0
end

################################################################################
# SphereCollider
################################################################################
mutable struct SphereCollider110{T} <: Collider{T}
    origin::Vector{T}
    radius::T
end

function SphereCollider(origin, radius)
    return SphereCollider110(origin, radius)
end

function inside(collider::SphereCollider110, p)
    origin = collider.origin
    return norm(p - origin) <= radius
end

################################################################################
# SoftCollider
################################################################################
mutable struct SoftCollider110{T,N} <: Collider{T}
    x::AbstractVector{T}
    q::Quaternion{T}
    mass::T
    inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    density_gradients::Vector{SVector{3,T}}
    weights::Vector{T} # contribution to collision force F = spring_constant * collision_weight
    weight_gradients::Vector{SVector{3,T}}
    nerf_object::Any
    mesh::GeometryBasics.Mesh
    options::ColliderOptions{T}
end

function SoftCollider(nerf_object, mesh; N=1000, density_scale=0.1, opts=ColliderOptions110(), T=Float64)
    x = szeros(T,3)
    q = Quaternion(1,0,0,0.0)
    mass, inertia, center_of_mass = inertia_properties(nerf_object, density_scale=density_scale)
    particles, densities, density_gradients = sample_soft(nerf_object, N)
    weights = densities ./ sum(densities) * mass
    weight_gradients = density_gradients ./ sum(densities) * mass
    return SoftCollider110{T,N}(x, q, mass, inertia, center_of_mass, particles,
        densities, density_gradients, weights, weight_gradients, nerf_object, mesh, opts)
end

function sample_soft(nerf_object, N::Int, T=Float64, min_density=1.0, max_density=105.0)
    particles = Vector{SVector{3,T}}(undef, N)
    densities = zeros(T,N)
    density_gradients = Vector{SVector{3,T}}(undef, N)

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
        density_gradients[i] = finite_difference_gradient(nerf_object, particles[i])
    end
    return particles, densities, density_gradients
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
function collision(collider::Collider{T}, soft_collider::SoftCollider110{T,N}) where {T,N}
    # soft.x = position of center of mass of soft in word frame
    # soft.q = orientation of body wrt world frame B -> qb -> W
    # P0 'origin' of the frame used by the nerf
    # Pcom center of mass of the body
    # pcom = position of Pcom in the nerf frame
    # particle = position of the particle in the nerf frame
    # p = position of the particle in the center of mass' frame attach to the body
    # p = x + qb * (particle - pcom)
    x = soft_collider.x
    q = soft_collider.q
    center_of_mass = soft_collider.center_of_mass
    ψ = Vector{T}()
    active_particles = []

    for i = 1:N
        particle = soft_collider.particles[i]
        p = x + Dojo.vector_rotate(particle - center_of_mass, q)
        if inside(collider, p)
            push!(ψ, soft_collider.weights[i])
            push!(active_particles, particle)
        end
    end

    contact_normals = [collision_normal(collider, soft_collider, p) for p in active_particles]
    num_active = length(ψ)
    return ψ, active_particles, contact_normals, num_active
end

function cross_collision(collider1::SoftCollider110{T,N1}, collider2::SoftCollider110{T,N2}) where {T,N1,N2}
    # returns
    # the particles of collider 1 in world frame
    # cross weights of collider1 and collider2
    # contact_normals for each particle
    # barycenter of contact of collider1 into collider2
    # contact normal from collider1 to collider2

    particles_w = zeros(Float32,N1,3)
    particles_2 = zeros(Float32,N1,3)
    for i = 1:N1
        pw = collider1.x + Dojo.vector_rotate(collider1.particles[i], collider1.q)
        p2 = Dojo.vector_rotate(pw - collider2.x, inv(collider2.q))
        particles_w[i,:] = pw
        particles_2[i,:] = p2
    end
    densities = py"density_query"(nerf_object, convert.(Float32, particles_2))
    weights = densities ./ sum(collider2.densities) * collider2.mass

    ψ = 0.0
    barycenter_w = zeros(3)
    contact_normal_w = zeros(3)
    for i = 1:N1
        weight_product = collider1.weights[i] * weights[i]
        ψ += weight_product
        barycenter_w += weight_product * particles_w[i,:]
        contact_normal_w += weight_product * Dojo.vector_rotate(-collider1.weight_gradients[i], collider1.q)
    end
    barycenter_w ./= (1e-20 + ψ)
    contact_normal_w /= (1e-20 + norm(contact_normal_w))
    return ψ, contact_normal_w, barycenter_w
end

function collision(collider1::SoftCollider110{T,N1}, collider2::SoftCollider110{T,N2}) where {T,N1,N2}
    # contact_normal = direction of force applied by collider1 on collider2
    ψ1, contact_normal1w, barycenter1w = cross_collision(collider1, collider2) # weigths of collider2 at particle1
    ψ2, contact_normal2w, barycenter2w = cross_collision(collider2, collider1) # weigths of collider1 at particle2
    ψ = ψ1 + ψ2
    contact_normal_w = (ψ1 * contact_normal1w + ψ2 * -contact_normal2w)
    contact_normal_w /= (1e-20 + norm(contact_normal_w))
    barycenter_w = (ψ1 * barycenter1w + ψ2 * barycenter2w) / (1e-20 + ψ)
    barycenter_2 = Dojo.vector_rotate(barycenter_w - collider2.x, inv(collider2.q))
    ∇ψ = zeros(3)
    return ψ, contact_normal_w, barycenter_2
end

function collision_normal(collider::HalfSpaceCollider110{T}, soft_collider::SoftCollider110{T,N}, particle) where {T,N}
    collider.normal
end

function collision_normal(collider::SphereCollider110{T}, soft_collider::SoftCollider110{T,N}, particle) where {T,N}
    normal = particle - collider.origin
    return normal ./ (1e-20 + norm(normal))
end
