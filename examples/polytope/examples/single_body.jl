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

# Nerf data
results_dir = joinpath(module_dir(), "results")
vertices = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["vertices"]
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
hull_vrep = removevredundancy(vrep(vertices), GLPK.Optimizer)
hull = polyhedron(hull_vrep)
hull_hrep = hrep(hull)

vis = Visualizer()
open(vis)

################################################################################
# Methods
################################################################################

function build_polytope!(poly::Polyhedron, vis::Visualizer; name::Symbol=:robot)
    mesh = Polyhedra.Mesh(poly)
    build_mesh!(mesh, vis; name=name)
end

function build_mesh!(mesh, vis::Visualizer; name::Symbol=:robot)
    setobject!(vis[name], mesh, MeshPhongMaterial(color=RGBA{Float32}(1, 1, 1, 0.7)))
    return nothing
end

function set_robot!(x, q, vis::Visualizer; name::Symbol=:robot)
    settransform!(vis[name], MeshCat.compose(
        MeshCat.Translation(x),
        MeshCat.LinearMap(q),)
        )
    return nothing
end

function object_intersection(x::AbstractVector, q::Quaternion, hrep::Polyhedra.Intersection)
    # x: position of the object expressed in the world frame
    # q: orientation of the object, q = object frame -> world frame
    z = SVector(0,0,1.0)
    zr = Dojo.vector_rotate(z, inv(q))
    b = z'*x
    floor_halfspace = HalfSpace(zr, -b)# in the frame of the object
    intersection_hrep = hull_hrep ∩ floor_halfspace
    intersection = polyhedron(intersection_hrep)
    vrep(intersection)
    return intersection
end

function contact_potential(x, q, hull::Polyhedron, nerf_object; N::Int=100, δx=0.002)
    intersection = object_intersection(x, q, hull_hrep)
    low_variance_sampling(intersection, nerf_object, N=N, δx=δx)
end

function low_variance_sampling(polytope, nerf_object; N::Int, δx=0.002)
    points = zeros(Float32, 6N, 3)
    vertices = polytope.vrep.points.points
    barycenter = mean(vertices)
    typical_length = mean(norm.(vertices .- [barycenter]))
    typical_surface = typical_length^2
    typical_volume = typical_length^3

    n = length(vertices)
    (n == 0) && return 0.0, zeros(3), [0,0,1]
    for i = 1:N
        p = (vertices[(i-1) % n + 1] + barycenter) / (2 + Int(floor(i/n)))
        points[(i-1)*6+1,:] = p + [-δx,0,0]
        points[(i-1)*6+2,:] = p + [+δx,0,0]
        points[(i-1)*6+3,:] = p + [0,-δx,0]
        points[(i-1)*6+4,:] = p + [0,+δx,0]
        points[(i-1)*6+5,:] = p + [0,0,-δx]
        points[(i-1)*6+6,:] = p + [0,0,+δx]
    end
    densities = py"density_query"(nerf_object, points)
    ϕ = typical_volume * mean(densities)

    ∇ϕ = zeros(3)
    for i = 1:N
        ∇ϕ += [
            densities[(i-1)*6+2] - densities[(i-1)*6+1],
            densities[(i-1)*6+4] - densities[(i-1)*6+3],
            densities[(i-1)*6+6] - densities[(i-1)*6+5]]
    end
    ∇ϕ *= typical_volume / (2δx) / N

    impact_normal = -mean(points, dims=1)[1,:]
    impact_normal /= norm(impact_normal)
    return ϕ, ∇ϕ, impact_normal
end

Dojo.@with_kw mutable struct SimulationOptions11{T}
    impact_damper::T=100.0
    impact_spring::T=30.0
    friction_coefficient::T=10.0
    friction_tanh::T=100.0
end

function integrator(y0; opts=SimulationOptions11())
    x = y0[1:3]
    v = y0[4:6]
    q = Quaternion(y0[7:10]..., false)
    ω = y0[11:13]

    ϕ, ∇ϕ, impact_normal = contact_potential(p0, q0, hull, nerf_object; N=100, δx=0.004)
    # ϕ = max(ϕ - 1, 0)
    impact_normal = [0,0,1]
    v_normal = v' * impact_normal * impact_normal
    v_tangential = v - v_normal

    # N = [g +
    #     1/m * impact_normal .* -opts.impact_damper * ϕ0 .* v0 +
    #     1/m * impact_normal * opts.impact_spring * ϕ0 +
    #     -1/m * opts.friction_coefficient * atan(opts.friction_tanh * norm(v0_tangential)) * v0_tangential/norm(v0_tangential) * ϕ0;
    #     zeros(3)]
    # y1 = (I - 0.5B) \ ((I + 0.5B) * y0) + A \ (N - (I - 0.5B) \ ((I + 0.5B) * N))
    F = g +
        1/m * impact_normal .* -opts.impact_damper * ϕ .* v +
        1/m * impact_normal * opts.impact_spring * ϕ +
        -1/m * opts.friction_coefficient * atan(opts.friction_tanh * norm(v_tangential)) * v_tangential/norm(v_tangential) * ϕ;
    τ = zeros(3)

    M = cat(m*Diagonal(ones(3)), J, dims=(1,2)) # mass matrix
    a = M \ ([F; τ] - [zeros(3); Dojo.skew(ω)* J * ω])
    x1 = x + h * v
    v1 = v + h * a[1:3]
    q1 = Dojo.next_orientation(q, SVector{3}(ω), h)
    ω1 = ω + h * a[4:6]
    y1 = [x1; v1; Dojo.vector(q1); ω1]
    return y1, ϕ, ∇ϕ
end

function exp_simulation(y0, N; opts=SimulationOptions11())
    Y = [y0]
    Φ = []
    ∇Φ = []
    for i = 1:N
        @show i
        y0, ϕ0, ∇ϕ0 = integrator(y0, opts=opts)
        push!(Y, y0)
        push!(Φ, ϕ0)
        push!(∇Φ, ∇ϕ0)
    end
    return Y, Φ, ∇Φ
end

set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.4))
set_background!(vis)


h = 0.01
m = 1.0
J = Diagonal(ones(3))
g = [0,0,-9.81]
A = [zeros(3,3) -1e-10I(3); -I(3) zeros(3,3)]
B = -A * h
H = 400

x0 = [0,0,1.0]
v0 = [0,0,0.0]
q0 = Quaternion(1,0,0,0.0,false)
ω0 = [0,0,0.0]
y0 = [x0; v0; Dojo.vector(q0); ω0]

opts = SimulationOptions11(
    impact_damper=100.0,
    impact_spring=30.0,
    friction_coefficient=5.0,
    friction_tanh=1000.0,
)
Y0, Φ0, ∇Φ0 = exp_simulation(y0, H, opts=opts)
plot([i*h for i=0:H], [p[6] for p in Y0], linewidth=3.0, xlabel="time", label="altitude")
plot([i*h for i=0:H], [p[5] for p in Y0], linewidth=3.0, xlabel="time", label="altitude")
plot!([i*h for i=1:H], Φ0, linewidth=3.0, xlabel="time", label="ϕ", color=:red)
# plot!([i*h for i=1:H], [p[3] for p in ∇Φ0], linewidth=3.0, xlabel="time", label="∇ϕ", color=:green)


setobject!(vis[:animated], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 1, 1, 1.0)))

animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    atframe(animation, i) do
        settransform!(vis[:animated], MeshCat.Translation(Y0[i][4:6]))
    end
end
setanimation!(vis, animation)




x0 = [0,0,0.0]
q0 = Quaternion(1,0,0,0.0)
contact_potential(x0, q0, hull, nerf_object, N=1000, δx=0.005)

ΔX = Vector(0.0001:0.0005:0.02)
plot(ΔX, [contact_potential(x0, q0, hull, nerf_object, N=1000, δx=δx)[2][1] for δx ∈ ΔX])

X = Vector(0.0:0.1:1.0)
plot(X, [contact_potential([0,0,x], q0, hull, nerf_object, N=1000, δx=0.004)[1] for x ∈ X])

X = -1.0:0.01:1.0
plot(X, 2 .* atan.(100*X)/π)

Dojo.convert_frames_to_video_and_gif("bunny_translation_slide")

build_polytope!(hull, vis)
build_mesh!(mesh, vis)
q = Quaternion(1,0.5,0,0.,true)
x = [0.1, 0.3, 0.9]
set_robot!(x, q, vis)

x = [1.2,0.4,0.4]
q = Quaternion(1,0.4,0,0.0,true)
intersection = object_intersection(hull_hrep, x, q)
build_polytope!(intersection, vis, name=:intersection)
set_robot!(x, q, vis, name=:intersection)
