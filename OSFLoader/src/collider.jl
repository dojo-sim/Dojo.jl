function inertia_properties(nerf_object; normalizer::DensityFieldNormalizer=DensityFieldNormalizer(),
        sampling_density::Int=30)
    load_python_script(mode=:cpu)

    xrange = range(-1.0, stop=1.0, length=sampling_density)
    yrange = range(-1.0, stop=1.0, length=sampling_density)
    zrange = range(-1.0, stop=1.0, length=sampling_density)
    particle_volume = (2/sampling_density)^3
    particles = grid_particles(xrange, yrange, zrange) # unnormalized particles

    densities = nerf_density(nerf_object, normalizer=normalizer, sampling_density=sampling_density)
    masses = normalizer.density_scale * particle_volume * densities

    mass = sum(masses)
    center_of_mass = zeros(3)
    for i = 1:sampling_density^3
        center_of_mass += masses[i] * particles[i,:]
    end
    center_of_mass /= mass
    inertia = zeros(3,3)
    for i = 1:sampling_density^3
        v = (particles[i,:] - center_of_mass)
        inertia += masses[i] * v * v'
    end
    return mass, inertia, center_of_mass
end

function finite_difference_gradient(nerf_object, particle;
        normalizer::DensityFieldNormalizer=DensityFieldNormalizer(), δ=0.01, n::Int=5)
    load_python_script(mode=:cpu)

    gradient = zeros(3)
    xrange = particle[1] .+ range(-δ, stop=δ, length=n)
    yrange = particle[2] .+ range(-δ, stop=δ, length=n)
    zrange = particle[3] .+ range(-δ, stop=δ, length=n)
    sample_particles = grid_particles(xrange, yrange, zrange)
    N = size(sample_particles)[1]
    for i = 1:N
        sample_particles[i,:] = particle_normalization(sample_particles[i,:]; normalizer=normalizer)
    end

    sample_gradients = py"density_gradient_query"(nerf_object, sample_particles)

    gradient = mean(sample_gradients, dims=1)[1,:]
    gradient /= normalizer.geometry_scaling
    R = rotationmatrix(normalizer.orientation_offset)
    gradient = R' * gradient
    return gradient
end

function sample_soft(nerf_object::Dict, num_particle::Int; T=Float64,
        normalizer::DensityFieldNormalizer=DensityFieldNormalizer(), particle_noise=0.005)

    density_fct = particles -> py"density_query"(nerf_object, particles)
    gradient_fct = particle -> finite_difference_gradient(nerf_object, particle)

    sample_soft(density_fct, gradient_fct, num_particle; T=T, normalizer=normalizer,
        particle_noise=particle_noise)
end

function sample_soft(density_fct::Function, gradient_fct::Function, num_particle::Int;
        T=Float64, normalizer::DensityFieldNormalizer=DensityFieldNormalizer(), particle_noise=0.005)
    load_python_script(mode=:cpu)
    density_low = normalizer.density_low
    density_high = normalizer.density_high

    particles = Vector{SVector{3,T}}(undef, num_particle)
    densities = zeros(T,num_particle)
    density_gradients = Vector{SVector{3,T}}(undef, num_particle)

    n = Int(floor((100num_particle)^(1/3)))
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    candidate_particles = grid_particles(xrange, yrange, zrange)
    candidate_particles += particle_noise * (rand(n^3, 3) .- 0.5)
    candidate_particles = convert.(Float32, candidate_particles)

    # normalizer
    normalized_particles = similar(candidate_particles)
    N = size(candidate_particles)[1]
    for i = 1:N
        normalized_particles[i,:] = particle_normalization(candidate_particles[i,:]; normalizer=normalizer)
    end

    normalized_densities = density_fct(normalized_particles)

    ind = findall(x -> density_low <= x <= density_high, normalized_densities)
    good_particles = candidate_particles[ind,:] # we take the unnormalized ones
    good_densities = normalized_densities[ind,:]
    num_good_particle = size(good_particles)[1]
    for i = 1:num_particle
        particles[i] = good_particles[Int(ceil(i*num_good_particle/num_particle)),:]
        densities[i] = good_densities[Int(ceil(i*num_good_particle/num_particle))]
        density_gradients[i] = gradient_fct(particles[i])
    end
    return particles, densities, density_gradients
end
