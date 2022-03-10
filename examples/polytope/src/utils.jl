function grid_points(xrange, yrange, zrange)
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

function slice_points(xrange, yrange, z)
    grid_points(xrange, yrange, z:z)
end

function grid_density(nerf_object, xrange, yrange, zrange)
    points = grid_points(xrange, yrange, zrange)
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

function naive_sample(poly, N::Int; rng::Int=0)
    Random.seed!(rng)
    points = zeros(Float32, N, 3)
    vertices = poly.vrep.points.points
    n = length(vertices)
    @assert n != 0
    Threads.@threads for i = 1:N
        θ = rand()
        points[i,:] = θ*vertices[rand(1:n)] + (1-θ)*vertices[rand(1:n)]
    end
    return points
end
