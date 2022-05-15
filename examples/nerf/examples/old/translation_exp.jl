################################################################################
# Translation only
################################################################################
using LinearAlgebra
using Plots
using FiniteDiff
using GLPK

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
nerf_object = py"generate_test_nerf"()

vis = Visualizer()
open(vis)

function get_intersection(hull_hrep, altitude)
    intersection_hrep = hull_hrep ∩ HalfSpace(SVector(0,0,1.0), -altitude)
    intersection = polyhedron(intersection_hrep)
    vrep(intersection)
    return intersection
end

function ray_sample(polytope, nerf_object, N::Int; rng::Int=0)
    Random.seed!(rng)
    points = zeros(Float32, N, 3)
    vertices = polytope.vrep.points.points
    barycenter = mean(vertices)
    typical_length = mean(norm.(vertices .- [barycenter]))
    typical_surface = typical_length^2
    typical_volume = typical_length^3

    n = length(vertices)
    Threads.@threads for i = 1:N
        θ = rand()
        points[i,:] = (vertices[(i-1) % n + 1] + barycenter) / (2 + Int(floor(i/n)))
    end
    ϕ = typical_volume * mean(py"density_query"(nerf_object, points))

    ∇ϕ = zeros(3)
    vertices_points = zeros(Float32, n, 3)
    for i = 1:n
        vertices_points[i,:] = vertices[i]
    end
    D = py"density_query"(nerf_object, vertices_points)
    for (i,v) in enumerate(vertices)
        d = v - barycenter
        d /= norm(d)
        ∇ϕ += typical_surface * D[i] * d
    end
    ∇ϕ = reshape(∇ϕ, 1,3)
    return ϕ, ∇ϕ
end

function interpenetration(x, nerf_object; N_samples::Int=100) where T
    # we trace N_samples rays between random pairs of vertices
    # we sample once along each ray, and we multpily by the length^3 of the ray
    # this means that larger intersections will have larger interpenetrations.
    altitude = x[3]
    intersection = get_intersection(hull_hrep, altitude)
    if length(intersection.vrep.points.points) > 0
        ϕ, ∇ϕ = ray_sample(intersection, nerf_object, N_samples)
    else
        points = []
        ϕ = 0.0
        ∇ϕ = zeros(1,3)
    end
    return ϕ, ∇ϕ
end

# Nerf data
results_dir = joinpath(module_dir(), "results")
vertices = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["vertices"]
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
hull_vrep = removevredundancy(vrep(vertices), GLPK.Optimizer)
hull = polyhedron(hull_vrep)
hull_hrep = hrep(hull)

interpenetration([0,0,0.201], nerf_object)
X = 0:0.0002:1.00
plot(X, [interpenetration([0,0,x], nerf_object, N_samples=1000)[1] for x in X])
plot(X, [interpenetration([0,0,x], nerf_object, N_samples=100)[3][1,3] for x in X])

δ = 0.01
X = 0:δ:1.00
D = [py"density_query"(nerf_object, [0 0 convert(Float32, x)])[1] for x in X]
G0 = (D[2:end] - D[1:end-1]) / δ
G = [py"density_gradient_query"(nerf_object, [0 0 convert(Float32, x)])[1,3] for x in X]

plot(X, D)
plot(X, G)
plot(X[1:end-1], G0)

py"density_gradient_query"(nerf_object, [0 0 convert(Float32, 0.1)])[1:3]
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100
a = 100


function soft_potential(x, nerf_object; N::Int=100)
    # we sample uniformly from the product of densities
        # normalized bunny density
        # soft 0-1 density for the halfspace
    altitude = x[3]
    intersection = get_intersection(hull_hrep, altitude)
    if length(intersection.vrep.points.points) > 0
        V = 0.0
        try
            V = Polyhedra.volume(intersection)
        catch e
        end
        points = naive_sample(intersection, N)
        ϕ = mean(py"density_query"(nerf_object, points)) / 100
        ∇ϕ = mean(py"density_gradient_query"(nerf_object, points), dims=1) / 100
        ϕ *= V
        ∇ϕ *= V
    else
        points = []
        ϕ = 0.0
        ∇ϕ = zeros(1,3)
    end
    return points, ϕ, ∇ϕ
end



function integrator(y0)
    v0 = y0[1:3]
    p0 = y0[4:6]
    points0, ϕ0, ∇ϕ0 = soft_potential(p0, nerf_object)
    N = [g + 1/m * [0,0,1] * 1000*ϕ0 + [0,0,1] * -100*∇ϕ0 * v0 ; zeros(3)]
    y1 = (I - 0.5B) \ ((I + 0.5B) * y0) + A \ (N - (I - 0.5B) \ ((I + 0.5B) * N))
    return y1
end

function exp_simulation(y0, N)
    Y = [y0]
    for i = 1:N
        @show i
        y0 = integrator(y0)
        push!(Y, y0)
    end
    return Y
end



# problem parameters
m = 1.0
g = -9.81
h = 0.01

p10 = [0,0,1.0]
p20 = [0,0,1.0]
ρ0 = 1.0
data0 = pack_data(p10, p20, ρ0)

p30 = [0,0,1.0]
γ0 = [0.0]
nγ = 1
vars0 = pack_vars(p30, γ0)

p0, ϕ0, ∇ϕ0 = soft_potential([0,0,0.3], nerf_object)
H = 40
residual(vars0, data0, ϕ0, ∇ϕ0)
soft_step(vars0, data0, ϕ0, ∇ϕ0)
P0, Γ0 = soft_simulate(vars0, data0, H)


setobject!(vis[:animated], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))

animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    P0[i]
    atframe(animation, i) do
        settransform!(vis[:animated], MeshCat.Translation(P0[i]))
    end
end
setanimation!(vis, animation)

plot([i*h for i=0:H+1], [p[3] for p in P0], linewidth=3.0, xlabel="time", ylabel="altitude")



intersection_hrep = hull_hrep ∩ HalfSpace(SVector(0,0,1.0),0.10)
intersection = polyhedron(intersection_hrep)
vrep(intersection)
intersection.vrep
alt_mesh = Polyhedra.Mesh(intersection)
setobject!(vis[:alt], alt_mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))






alts = 0:0.002:1
plot(alts, [soft_potential([0,0,alt], nerf_object, N=1000)[2] for alt in alts])
plot(alts, [soft_potential([0,0,alt], nerf_object, N=100)[3][1] for alt in alts])
plot(alts, [soft_potential([0,0,alt], nerf_object, N=100)[3][3] for alt in alts])

soft_potential([0,0,1.0], nerf_object, N=10)[2]
include("../src/utils.jl")

Dojo.convert_frames_to_video_and_gif("bunny_unstable_explicit_scheme")
using Dojo

_, _, ∇ϕ0 = soft_potential(p0, nerf_object)


h = 0.02
m = 1.0
g = [0,0,-9.81]
A = [zeros(3,3) -1e-10I(3); -I(3) zeros(3,3)]
B = -A * h
H = 400

v0 = [0,0,0.0]
p0 = [0,0,1.0]
y0 = [v0; p0]

Y0 = exp_simulation(y0, H)
plot([i*h for i=0:H], [p[6] for p in Y0], linewidth=3.0, xlabel="time", ylabel="altitude")

setobject!(vis[:animated], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(0, 0, 1, 0.5)))

animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    atframe(animation, i) do
        settransform!(vis[:animated], MeshCat.Translation(Y0[i][4:6]))
    end
end
setanimation!(vis, animation)
