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
@pyinclude(joinpath(OSF_PATH, "extract_density_julia.py"))
nerf_object = py"generate_test_nerf"()

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

@time collision(halfspace, soft)

################################################################################
# Methods
################################################################################
Dojo.@with_kw mutable struct SimulationOptions14{T}
    impact_damper::Vector{T}=[1,1,400.0]
    impact_spring::T=100.0
    friction_coefficient::T=10.0
    friction_tanh::T=100.0
    angular_damping::T=1.0
    angular_friction::T=1.0
end

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

    y1 = [x1; v1; Dojo.vector(q1); ω1]
    return y1, ϕ
end

function inegrator_line_search(a, res, Δ, residual)
    α = 1.0
    for i = 1:10
        (norm(residual(a + α * Δ), Inf) <= norm(res, Inf)) && break
        α *= 0.5
    end
    return α
end

function implicit_simulation(y0, N; opts=SimulationOptions14())
    y = [y0]
    ϕ = []
    for i = 1:N
        @show i
        y0, ϕ0 = implicit_integrator(y0, opts=opts)
        push!(y, y0)
        push!(ϕ, ϕ0)
    end
    return y, ϕ
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
y0 = [x0; v0; Dojo.vector(q0); ω0]

opts = SimulationOptions14(
    impact_damper=[1,1,500.0],
    impact_spring=100.0,
    friction_coefficient=0.2,
    friction_tanh=1000.0,
    angular_damping=0.5,
    angular_friction=0.04,
)
@time Y0, Φ0 = implicit_simulation(y0, H, opts=opts)
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
    end
end
setanimation!(vis, animation)
# render_static(vis)
# open(joinpath(results_dir, "bunny_hull.html"), "w") do file
#     write(file, static_html(vis))
# end
# Dojo.convert_frames_to_video_and_gif("bunny_rolling_friction_minus")
