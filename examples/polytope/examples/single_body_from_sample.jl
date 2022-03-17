################################################################################
# Translation only
################################################################################
using LinearAlgebra
using Plots
using FiniteDiff
using GLPK
using Dojo
using PyCall

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
nerf_object = py"generate_test_nerf"()

# Nerf data
results_dir = joinpath(example_dir(), "results")
mesh_vertices = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["vertices"]
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
tight_mesh = jldopen(joinpath(results_dir, "bunny_tight_mesh.jld2"))["mesh"]
hull_vrep = removevredundancy(vrep(mesh_vertices), GLPK.Optimizer)
hull = polyhedron(hull_vrep)
hull_hrep = hrep(hull)

vis = Visualizer()
open(vis)

################################################################################
# Offline
################################################################################

abstract type Collider{T} end

mutable struct HalfSpaceCollider25{T} <: Collider{T}
    altitude::T
end

mutable struct SoftCollider25{T,N} <: Collider{T}
    x::AbstractVector{T}
    q::Quaternion{T}
    mass::T
    inertia::AbstractMatrix{T}
    center_of_mass::AbstractVector{T}
    particles::Vector{SVector{3,T}}
    densities::Vector{T}
    gradients::Vector{SVector{3,T}}
    nerf_object::Any
    mesh::GeometryBasics.Mesh
end

function HalfSpaceCollider(altitude)
    return HalfSpaceCollider25(altitude)
end

function SoftCollider(nerf_object, mesh; N=1000, T=Float64)
    x = szeros(T,3)
    q = Quaternion(1,0,0,0.0)
    mass, inertia, center_of_mass = inertia_properties(nerf_object)
    particles, densities, gradients = sample_soft(nerf_object, N)
    return SoftCollider25{T,N}(x, q, mass, inertia, center_of_mass, particles, densities, gradients, nerf_object, mesh)
end

function inertia_properties(nerf_object; n=30, scale=1e-1)
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    particle_volume = (2/n)^3
    particles = grid_points(xrange, yrange, zrange)
    masses = particle_volume * scale .* py"density_query"(nerf_object, particles)

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

function sample_soft(nerf_object, N::Int, T=Float64, min_density=1.0, max_density=105.0)
    particles = Vector{SVector{3,T}}(undef, N)
    densities = zeros(T,N)
    gradients = Vector{SVector{3,T}}(undef, N)

    n = Int(floor((100N)^(1/3)))
    xrange = range(-1.0, stop=1.0, length=n)
    yrange = range(-1.0, stop=1.0, length=n)
    zrange = range(-1.0, stop=1.0, length=n)
    candidate_particles = grid_points(xrange, yrange, zrange)
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
    sample_particles = grid_points(xrange, yrange, zrange)
    sample_gradients = py"density_gradient_query"(nerf_object, sample_particles)
    gradient = mean(sample_gradients, dims=1)[1,:]
    return gradient
end

function build_collider!(collider::SoftCollider25{N,T}, vis::Visualizer;
        name::Symbol=:collider, visualize_particle::Bool=false,
        mesh_color=RGBA(1,1,1,0.3),
        particle_color=RGBA(0.8,0.8,0.8,1.0),
        gradient_color=RGBA(0.9,0.2,0.2,1.0),
        com_color=RGBA(0,0,0,1.0),
        ) where {T,N}
    setobject!(vis[name][:mesh], collider.mesh, MeshPhongMaterial(color=mesh_color))
    setobject!(vis[name][:center_of_mass], HyperSphere(Point(collider.center_of_mass...), 0.03), MeshPhongMaterial(color=com_color))
    if visualize_particle
        for i = 1:N
            particle = collider.particles[i]
            gradient = collider.gradients[i]
            setobject!(vis[name][:particles]["$i"], HyperSphere(Point(particle...), 0.005), MeshPhongMaterial(color=particle_color))
            arrow_vis = ArrowVisualizer(vis[name][:gradients]["$i"])
            setobject!(arrow_vis, MeshPhongMaterial(color=gradient_color))
            settransform!(arrow_vis, Point(particle...), Vec(-2e-5gradient...), shaft_radius=0.0025)
        end
    end
    return nothing
end

function collision(halfspace::HalfSpaceCollider25{T}, collider::SoftCollider25{T,N}) where {T,N}
    ϕ = 0.0
    ∇ϕ = zeros(3)
    contact_normal = [0,0,1.0]
    barycenter = zeros(3)

    active = 0
    for i = 1:N
        particle = collider.particles[i]
        p = collider.x + Dojo.vector_rotate(particle, collider.q)
        if p[3] < halfspace.altitude
            active +=1
            density = collider.densities[i]
            gradient = collider.gradients[i]
            normal = gradient / norm(gradient)
            ϕ += density
            ∇ϕ += density * gradient
            barycenter += density * particle
        end
    end
    if active > 0
        ∇ϕ /= ϕ * N
        barycenter /= ϕ
        # contact_normal = ∇ϕ / norm(∇ϕ)
        ϕ /= N
    end
    return ϕ, ∇ϕ, contact_normal, barycenter
end


altitude = 0.0
halfspace = HalfSpaceCollider(altitude)
sample_soft(nerf_object, N)
# soft = SoftCollider(nerf_object, mesh, N=5000)

soft.x = [0,0,0.3]
@time collision(halfspace, soft)

X = -1.0:0.01:2.0
phi = zeros(length(X))
∇phi = zeros(length(X),3)
Nphi = zeros(length(X),3)
for (i,x) in enumerate(X)
    soft.x = [0,0,x]
    phi[i], ∇phi[i,:], Nphi[i,:] = collision(halfspace, soft)
end
plot(X, phi)
plot!(X, ∇phi[:,1])
plot!(X, ∇phi[:,2])
plot!(X, ∇phi[:,3])
plot!(X, Nphi[:,1])
plot!(X, Nphi[:,2])
plot!(X, Nphi[:,3])



# vis = Visualizer()
# open(vis)
# render_static(vis)
# open(joinpath(results_dir, "bunny_gradients.html"), "w") do file
#     write(file, static_html(vis))
# end

set_light!(vis)
set_background!(vis)
set_floor!(vis, alt=-0.5)
build_collider!(soft, vis, visualize_particle=false)




################################################################################
# Methods
################################################################################
Dojo.@with_kw mutable struct SimulationOptions14{T}
    impact_damper::Vector{T}=[1,1,100.0]
    impact_spring::T=30.0
    friction_coefficient::T=10.0
    friction_tanh::T=100.0
    angular_damping::T=1.0
    angular_friction::T=1.0
end

# function integrator(y0; opts=SimulationOptions14())
#     x = y0[1:3]
#     v = y0[4:6]
#     q = Quaternion(y0[7:10]..., false)
#     ω = y0[11:13]
#
#     soft.x = x
#     soft.q = q
#     ϕ, ∇ϕ, impact_normal, barycenter = collision(halfspace, soft)
#
#     barycenter_w = Dojo.vector_rotate(barycenter, q)
#     vc = v + Dojo.vector_rotate(Dojo.skew(-barycenter) * ω, q)
#     vc_normal = vc' * impact_normal * impact_normal
#     vc_tangential = vc - vc_normal
#
#     F_impact = -opts.impact_damper * ϕ .* vc
#     F_impact += impact_normal * opts.impact_spring * ϕ
#     F_friction = -opts.friction_coefficient * atan(opts.friction_tanh * norm(vc_tangential)) * vc_tangential/(1e-3 + norm(vc_tangential)) * norm(F_impact)
#     F_contact = F_impact + F_friction
#     F = m*g + F_contact
#
#     τ_w = Dojo.skew(barycenter_w) * F_contact
#     τ = Dojo.vector_rotate(τ_w, inv(q))
#     τ -= opts.angular_damping * ω  * ϕ
#     τ -= opts.angular_friction * atan(opts.friction_tanh * norm(ω)) * ω/(1e-3 + norm(ω)) * norm(F_impact)
#
#     M = cat(m*Diagonal(ones(3)), J, dims=(1,2)) # mass matrix
#     a = M \ ([F; τ] - [zeros(3); Dojo.skew(ω)* J * ω])
#
#     v1 = v + h * a[1:3]
#     x1 = x + h * (v+v1)/2
#     ω1 = ω + h * a[4:6]
#     q1 = Dojo.next_orientation(q, SVector{3}((ω+ω1)/2), h)
#     y1 = [x1; v1; Dojo.vector(q1); ω1]
#     return y1, ϕ, ∇ϕ, τ
# end
#

function implicit_integrator(y0; opts=SimulationOptions14())
    x = y0[1:3]
    v = y0[4:6]
    q = Quaternion(y0[7:10]..., false)
    ω = y0[11:13]

    soft.x = x
    soft.q = q
    ϕ, ∇ϕ, impact_normal, barycenter = collision(halfspace, soft)

    function residual(a)
        v1 = v + h * a[1:3]
        ω1 = ω + h * a[4:6]
        x1 = x + h * (v+v1)/2
        q1 = Dojo.next_orientation(q, SVector{3}((ω+ω1)/2), h)

        barycenter_w = Dojo.vector_rotate(barycenter, q) # TODO could be q1
        vc = v1 + Dojo.vector_rotate(Dojo.skew(-barycenter) * ω1, q) # TODO could be q1
        vc_normal = vc' * impact_normal * impact_normal
        vc_tangential = vc - vc_normal

        F_impact = -opts.impact_damper * ϕ .* vc
        F_impact += impact_normal * opts.impact_spring * ϕ
        F_friction = -opts.friction_coefficient * atan(opts.friction_tanh * norm(vc_tangential)) * vc_tangential/(1e-3 + norm(vc_tangential)) * norm(F_impact)
        F_contact = F_impact + F_friction
        F = m*g + F_contact

        τ_w = Dojo.skew(barycenter_w) * F_contact
        τ = Dojo.vector_rotate(τ_w, inv(q)) # TODO could be q1
        τ -= opts.angular_damping * ω1 * ϕ
        τ -= opts.angular_friction * atan(opts.friction_tanh * norm(ω1)) * ω1/(1e-3 + norm(ω1)) * norm(F_impact)

        M = cat(m*Diagonal(ones(3)), J, dims=(1,2)) # mass matrix
        ā = M \ ([F; τ] - [zeros(3); Dojo.skew(ω1) * J * ω1])
        return ā - a
    end

    function jacobian_residual(a)
        jac = FiniteDiff.finite_difference_jacobian(a -> residual(a), a)
        return jac
    end

    function solve_integrator()
        a = zeros(6)
        for i = 1:10
            res = residual(a)
            @show i
            @show norm(res, Inf)
            (norm(res, Inf) < 1e-4) && break
            jac = jacobian_residual(a)
            Δ = - jac \ res
            α = inegrator_line_search(a, res, Δ, residual)
            a += α * Δ
        end
        return a
    end
    a = solve_integrator()
    v1 = v + h * a[1:3]
    ω1 = ω + h * a[4:6]
    x1 = x + h * (v+v1)/2
    q1 = Dojo.next_orientation(q, SVector{3}((ω+ω1)/2), h)

    y1 = [x1; v1; Dojo.vector(q1); ω1]
    return y1, ϕ, ∇ϕ#, τ
end

function inegrator_line_search(a, res, Δ, residual)
    α = 1.0
    for i = 1:10
        (norm(residual(a + α * Δ), Inf) <= norm(res, Inf)) && break
        α *= 0.5
    end
    return α
end

implicit_integrator(y0; opts=SimulationOptions14())


function implicit_simulation(y0, N; opts=SimulationOptions14())
    y = [y0]
    ϕ = []
    ∇ϕ = []
    for i = 1:N
        @show i
        y0, ϕ0, ∇ϕ0= implicit_integrator(y0, opts=opts)
        push!(y, y0)
        push!(ϕ, ϕ0)
        push!(∇ϕ, ∇ϕ0)
    end
    return y, ϕ, ∇ϕ
end

# function explicit_simulation(y0, N; opts=SimulationOptions14())
#     y = [y0]
#     ϕ = []
#     ∇ϕ = []
#     τ = []
#     for i = 1:N
#         @show i
#         y0, ϕ0, ∇ϕ0, τ0 = integrator(y0, opts=opts)
#         push!(y, y0)
#         push!(ϕ, ϕ0)
#         push!(∇ϕ, ∇ϕ0)
#         push!(τ, τ0)
#     end
#     return y, ϕ, ∇ϕ, τ
# end

set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5)
set_background!(vis)

h = 0.02
m = soft.mass
J = soft.inertia
g = [0,0,-9.81]
A = [zeros(3,3) -1e-10I(3); -I(3) zeros(3,3)]
B = -A * h
H = Int(floor(3.0/h))

x0 = [0,0,1.0]
v0 = [0,2,3.0]
q0 = rand(4)
q0 ./= norm(q0)
q0 = Quaternion(q0...)
ω0 = [3.5,0,0.0]
y0 = [x0; v0; Dojo.vector(q0); ω0]

opts = SimulationOptions14(
    impact_damper=[1,1,400.0],
    impact_spring=100.0,
    friction_coefficient=0.1,
    friction_tanh=1000.0,
    angular_damping=0.5,
    angular_friction=0.04,
)
Y0, Φ0, ∇Φ0 = implicit_simulation(y0, H, opts=opts)
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



build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(0.6,0.8,0.9,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
# build_polytope!(hull, vis, name=:hull, color=RGBA(1,0.3,0.3,0.3))
animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    atframe(animation, i) do
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:mesh)
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:tight_mesh)
        set_robot!(Y0[i][1:3], Quaternion(Y0[i][7:10]...), vis, name=:hull)
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
# Dojo.convert_frames_to_video_and_gif("bunny_hull_mesh_tight_mesh")
