function grid_particles(xrange, yrange, zrange)
    nx = length(xrange)
    ny = length(yrange)
    nz = length(zrange)
    xv = Vector(xrange)
    yv = Vector(yrange)
    zv = Vector(zrange)
    points = zeros(Float32, nx*ny*nz, 3)
    for i = 1:nx
        for j = 1:ny
            for k = 1:nz
                points[1+(i-1)+(j-1)*nx+(k-1)*ny*nx,:] .= [xv[i], yv[j], zv[k]]
            end
        end
    end
    return points
end

function slice_particles(xrange, yrange, z)
    grid_particles(xrange, yrange, z:z)
end

function grid_density(nerf_object, xrange, yrange, zrange)
    points = grid_particles(xrange, yrange, zrange)
    density = py"density_query"(nerf_object, points)
    nx = length(xrange)
    ny = length(yrange)
    nz = length(zrange)
    density = reshape(density, nx, ny, nz)
    return density
end

function slice_density(nerf_object, xrange, yrange, z)
    grid_density(nerf_object, xrange, yrange, z:z)[:,:,1]
end

function inertia_properties(nerf_object; n=30, density_scale=1e-1)
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    particle_volume = (2/n)^3
    particles = grid_particles(xrange, yrange, zrange)
    masses = density_scale * particle_volume * py"density_query"(nerf_object, particles)

    mass = sum(masses)
    center_of_mass = zeros(3)
    for i = 1:n^3
        center_of_mass += masses[i] * particles[i,:]
    end
    center_of_mass /= mass
    inertia = zeros(3,3)
    for i = 1:n^3
        v = (particles[i,:] - center_of_mass)
        inertia += masses[i] * v * v'
    end
    return mass, inertia, center_of_mass
end


function sample_soft(nerf_object, N::Int; T=Float64, min_density=1.0, max_density=105.0,
        particle_noise=0.005)
    particles = Vector{SVector{3,T}}(undef, N)
    densities = zeros(T,N)
    density_gradients = Vector{SVector{3,T}}(undef, N)

    n = Int(floor((100N)^(1/3)))
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    candidate_particles = grid_particles(xrange, yrange, zrange)
    candidate_particles += particle_noise * (rand(n^3, 3) .- 0.5)
    candidate_particles = convert.(Float32, candidate_particles)
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
