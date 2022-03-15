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
tight_mesh = jldopen(joinpath(results_dir, "bunny_tight_mesh.jld2"))["mesh"]
hull_vrep = removevredundancy(vrep(vertices), GLPK.Optimizer)
hull = polyhedron(hull_vrep)
hull_hrep = hrep(hull)

vis = Visualizer()
open(vis)

################################################################################
# Methods
################################################################################

function build_polytope!(poly::Polyhedron, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    mesh = Polyhedra.Mesh(poly)
    build_mesh!(mesh, vis; name=name, color=color)
end

function build_mesh!(mesh, vis::Visualizer; name::Symbol=:robot, color=RGBA(1,1,1,1.0))
    setobject!(vis[name], mesh, MeshPhongMaterial(color=color))
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
    (n == 0) && return 0.0, zeros(3), [0,0,1], [0,0,0]
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
    return ϕ, ∇ϕ, impact_normal, barycenter
end

Dojo.@with_kw mutable struct SimulationOptions14{T}
    impact_damper::Vector{T}=[1,1,100.0]
    impact_spring::T=30.0
    friction_coefficient::T=10.0
    friction_tanh::T=100.0
    angular_damping::T=1.0
    angular_friction::T=1.0
end

function integrator(y0; opts=SimulationOptions14())
    x = y0[1:3]
    v = y0[4:6]
    q = Quaternion(y0[7:10]..., false)
    ω = y0[11:13]

    ϕ, ∇ϕ, impact_normal, barycenter = contact_potential(x, q, hull, nerf_object; N=100, δx=0.004)
    ϕimp = max(ϕ-0.5,0)
    impact_normal = [0,0,1]
    barycenter_w = Dojo.vector_rotate(barycenter, q)
    vc = v + Dojo.vector_rotate(Dojo.skew(-barycenter) * ω, q)
    vc_normal = vc' * impact_normal * impact_normal
    vc_tangential = vc - vc_normal

    # F_impact = impact_normal .* -opts.impact_damper * ϕ .* vc
    F_impact = -opts.impact_damper * ϕ .* vc
    F_impact += impact_normal * opts.impact_spring * ϕimp^2/(3+ϕimp)
    F_impact += -opts.friction_coefficient * atan(opts.friction_tanh * norm(vc_tangential)) * vc_tangential/(1e-3 + norm(vc_tangential)) * ϕ
    F = m*g + F_impact

    τ_w = Dojo.skew(barycenter_w) * F_impact
    τ = Dojo.vector_rotate(τ_w, inv(q))
    τ -= opts.angular_damping * ω  * ϕ
    τ -= opts.angular_friction * atan(opts.friction_tanh * norm(ω)) * ω/(1e-3 + norm(ω)) * ϕ

    M = cat(m*Diagonal(ones(3)), J, dims=(1,2)) # mass matrix
    a = M \ ([F; τ] - [zeros(3); Dojo.skew(ω)* J * ω])

    v1 = v + h * a[1:3]
    x1 = x + h * (v+v1)/2
    ω1 = ω + h * a[4:6]
    q1 = Dojo.next_orientation(q, SVector{3}((ω+ω1)/2), h)
    y1 = [x1; v1; Dojo.vector(q1); ω1]
    return y1, ϕ, ∇ϕ, τ
end

function exp_simulation(y0, N; opts=SimulationOptions14())
    y = [y0]
    ϕ = []
    ∇ϕ = []
    τ = []
    for i = 1:N
        @show i
        y0, ϕ0, ∇ϕ0, τ0 = integrator(y0, opts=opts)
        push!(y, y0)
        push!(ϕ, ϕ0)
        push!(∇ϕ, ∇ϕ0)
        push!(τ, τ0)
    end
    return y, ϕ, ∇ϕ, τ
end

set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5)
set_background!(vis)


h = 0.005
m = 1.0
J = 1Diagonal(ones(3))
g = [0,0,-9.81]
A = [zeros(3,3) -1e-10I(3); -I(3) zeros(3,3)]
B = -A * h
H = 1000

x0 = [0,0,1.0]
v0 = [0,2,3.0]
q0 = rand(4)
q0 ./= norm(q0)
q0 = Quaternion(q0...)
ω0 = [1.5,0,0.0]
y0 = [x0; v0; Dojo.vector(q0); ω0]

opts = SimulationOptions14(
    impact_damper=[1,1,60.0],
    impact_spring=20.0,
    friction_coefficient=1.0,
    friction_tanh=1000.0,
    angular_damping=1.0,
    angular_friction=0.2,
)
Y0, Φ0, ∇Φ0 = exp_simulation(y0, H, opts=opts)
plt = plot(layout=(3,1))
plot!(plt[1,1], [i*h for i=0:H], [y[1] for y in Y0], linewidth=3.0, xlabel="time", label="x", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [y[2] for y in Y0], linewidth=3.0, xlabel="time", label="y", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [y[3] for y in Y0], linewidth=3.0, xlabel="time", label="z", color=:blue)
plot!(plt[2,1], [i*h for i=0:H], [y[7] for y in Y0], linewidth=3.0, xlabel="time", label="q1", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [y[8] for y in Y0], linewidth=3.0, xlabel="time", label="q2", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [y[9] for y in Y0], linewidth=3.0, xlabel="time", label="q3", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [y[10] for y in Y0], linewidth=3.0, xlabel="time", label="q4", color=:green)
plot!(plt[3,1], [i*h for i=1:H], Φ0, linewidth=3.0, xlabel="time", label="ϕ", color=:red)
# plot!([i*h for i=1:H], [p[3] for p in ∇Φ0], linewidth=3.0, xlabel="time", label="∇ϕ", color=:green)

setobject!(vis[:animated], mesh,
    MeshPhongMaterial(color=RGBA{Float32}(1, 1, 1, 1.0)))

build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.4))
build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(0.4,0.0,0.8,0.6))
build_polytope!(hull, vis, name=:hull, color=RGBA(1,0.3,0.3,0.6))
animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    atframe(animation, i) do
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:mesh)
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:tight_mesh)
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:hull)
    end
end
setanimation!(vis, animation)







x0 = [0,0,0.0]
q0 = Quaternion(1,0,0,0.0)
contact_potential(x0, q0, hull, nerf_object, N=1000, δx=0.005)

ΔX = Vector(0.0001:0.0005:0.02)
plot(ΔX, [contact_potential(x0, q0, hull, nerf_object, N=1000, δx=δx)[2][1] for δx ∈ ΔX])

X = Vector(0.0:0.0002:1.0)
plot(X, [contact_potential([0,0,x], q0, hull, nerf_object, N=100, δx=0.004)[1] for x ∈ X])

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
