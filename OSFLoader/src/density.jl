function get_nerf_object(;osf_path=OSFLoader.OSF_PATH,
        config_folder=OSFLoader.OSF_CONFIG_FOLDER, filename="bunny_trans")
    config_path = joinpath(osf_path, config_folder, filename, filename *".txt")
    load_python_script(mode=:cpu)
    nerf_object = OSFLoader.py"generate_test_nerf"(config_path)
    return nerf_object
end

function density_query(nerf_object, particles)
    load_python_script(mode=:cpu)
    densities = py"density_query"(nerf_object, convert.(Float32, particles))
    return densities
end

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

function grid_density(nerf_object, xrange, yrange, zrange; normalizer=DensityFieldNormalizer())
    load_python_script(mode=:cpu)
    points = grid_particles(xrange, yrange, zrange)
    N = size(points)[1]
    for i = 1:N
        points[i,:] = particle_normalization(points[i,:]; normalizer=normalizer)
    end
    density = py"density_query"(nerf_object, points)
    nx = length(xrange)
    ny = length(yrange)
    nz = length(zrange)
    density = reshape(density, nx, ny, nz)
    return density
end

function slice_density(nerf_object, xrange, yrange, z; normalizer=DensityFieldNormalizer())
    grid_density(nerf_object, xrange, yrange, z:z; normalizer=normalizer)[:,:,1]
end

function nerf_density(nerf_object;
    normalizer::DensityFieldNormalizer=DensityFieldNormalizer(),
    sampling_density::Int=20)
    xrange = range(-1.0, stop=1.0, length=sampling_density)
    yrange = range(-1.0, stop=1.0, length=sampling_density)
    zrange = range(-1.0, stop=1.0, length=sampling_density)

    D = grid_density(nerf_object, xrange, yrange, zrange,
            normalizer=normalizer)
    return D
end
