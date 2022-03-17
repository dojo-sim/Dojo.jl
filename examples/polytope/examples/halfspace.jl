################################################################################
# Translation only
################################################################################
using LinearAlgebra
using Plots
using FiniteDiff
using Dojo
using PyCall

################################################################################
# build NeRf
################################################################################
# @pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
# nerf_object = py"generate_test_nerf"()

# Nerf data
results_dir = joinpath(example_dir(), "results")
mesh = jldopen(joinpath(results_dir, "bunny_mesh.jld2"))["mesh"]
tight_mesh = jldopen(joinpath(results_dir, "bunny_tight_mesh.jld2"))["mesh"]

vis = Visualizer()
open(vis)

################################################################################
# Offline
################################################################################
halfspace_origin = [0,0,0.0]
normal = [0.0,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)
# soft = SoftCollider(nerf_object, mesh, N=5000)
soft = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]
@time collision(halfspace, soft)

################################################################################
# Methods
################################################################################

function implicit_integrator(z)
    opts = soft.options
    x = z[1:3]
    v = z[4:6]
    q = Quaternion(z[7:10]..., false)
    ω = z[11:13]

    soft.x = x
    soft.q = q
    ϕ, ∇ϕ, impact_normal, barycenter = collision(halfspace, soft)

    coulomb_direction(v) = - atan(opts.coulomb_smoothing * norm(v)) * v/(opts.coulomb_regularizer + norm(v))

    function residual(a)
        v1 = v + h * a[1:3]
        ω1 = ω + h * a[4:6]
        x1 = x + h * (v+v1)/2
        q1 = Dojo.next_orientation(q, SVector{3}((ω+ω1)/2), h)

        barycenter_w = Dojo.vector_rotate(barycenter, q) # TODO could be q1
        vc = v1 + Dojo.vector_rotate(Dojo.skew(-barycenter) * ω1, q) # TODO could be q1
        vc_normal = vc' * impact_normal * impact_normal
        vc_tangential = vc - vc_normal

        F_impact = -opts.impact_damper * ϕ * vc
        F_impact += impact_normal * opts.impact_spring * ϕ
        F_friction = opts.sliding_friction * coulomb_direction(vc_tangential) * norm(F_impact)
        F_contact = F_impact + F_friction
        F = m*g + F_contact

        τ_w = Dojo.skew(barycenter_w) * F_contact
        τ = Dojo.vector_rotate(τ_w, inv(q)) # TODO could be q1
        τ += -opts.rolling_drag * ϕ * ω1
        τ += opts.rolling_friction * coulomb_direction(ω1) * norm(F_impact)

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
            (norm(res, Inf) < 1e-6) && break
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

    z1 = [x1; v1; Dojo.vector(q1); ω1]
    return z1, ϕ
end

function inegrator_line_search(a, res, Δ, residual)
    α = 1.0
    for i = 1:10
        (norm(residual(a + α * Δ), Inf) <= norm(res, Inf)) && break
        α *= 0.5
    end
    return α
end

function implicit_simulation(z0, N)
    z = [z0]
    ϕ = []
    for i = 1:N
        @show i
        z0, ϕ0 = implicit_integrator(z0)
        push!(z, z0)
        push!(ϕ, ϕ0)
    end
    return z, ϕ
end


set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5, alt=0.0)
set_background!(vis)

h = 0.01
H = Int(floor(4.0/h))
m = soft.mass
J = soft.inertia
g = [0,0,-9.81]

x0 = [0,0,1.0]
v0 = [0,2,3.0]
q0 = rand(4)
q0 ./= norm(q0)
q0 = Quaternion(q0...)
ω0 = [3.5,0,0.0]
z0 = [x0; v0; Dojo.vector(q0); ω0]

soft.options = ColliderOptions18(
    impact_damper=300.0,
    impact_spring=100.0,
    sliding_friction=0.2,
    coulomb_smoothing=1000.0,
    rolling_drag=0.5,
    rolling_friction=0.04,
)
@time Z0, Φ0 = implicit_simulation(z0, H)
plt = plot(layout=(3,1), xlabel="time")
plot!(plt[1,1], [i*h for i=0:H], [z[1] for z in Z0], linewidth=3.0, label="x", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[2] for z in Z0], linewidth=3.0, label="y", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[3] for z in Z0], linewidth=3.0, label="z", color=:blue)
plot!(plt[2,1], [i*h for i=0:H], [z[7] for z in Z0], linewidth=3.0, label="q1", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[8] for z in Z0], linewidth=3.0, label="q2", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[9] for z in Z0], linewidth=3.0, label="q3", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[10] for z in Z0], linewidth=3.0, label="q4", color=:green)
plot!(plt[3,1], [i*h for i=1:H], Φ0, linewidth=3.0, label="ϕ", color=:red)



build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(0.6,0.8,0.9,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
# build_polytope!(hull, vis, name=:hull, color=RGBA(1,0.3,0.3,0.3))
animation = MeshCat.Animation(Int(floor(1/h)))
for i = 1:H
    atframe(animation, i) do
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:mesh)
        set_robot!(Z0[i][1:3], Quaternion(Z0[i][7:10]...), vis, name=:tight_mesh)
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
# Dojo.convert_frames_to_video_and_gif("bunny_rolling_friction_minus")
