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
# Methods
################################################################################
function implicit_integrator(halfspace, soft, timestep, z)
    opts = soft.options
	x2, v15, q2, ϕ15 = unpack_maximal_state(z, 1)
	x1 = Dojo.next_position(x2, -v15, timestep)
	q1 = Dojo.next_orientation(q2, -ϕ15, timestep)

    soft.x = x2
    soft.q = q2
    ψ, ∇ψ, impact_normal, barycenter = collision(halfspace, soft)
    coulomb_direction(v) = - atan(opts.coulomb_smoothing * norm(v)) * v/(opts.coulomb_regularizer + norm(v))

    function residual(vϕ)
        v25 = vϕ[SUnitRange(1,3)]
        ϕ25 = vϕ[SUnitRange(4,6)]
        x3 = Dojo.next_position(x2, v25, timestep)
        q3 = Dojo.next_orientation(q2, ϕ25, timestep)

        # dynamics
        D1x = - 1.0 / timestep * mass * (x2 - x1) - 0.5 * timestep * mass * gravity
        D2x =   1.0 / timestep * mass * (x3 - x2) - 0.5 * timestep * mass * gravity
        D1q = -2.0 / timestep * LVᵀmat(q2)' * Lmat(q1) * Vᵀmat() * inertia * Vmat() * Lmat(q1)' * vector(q2)
        D2q = -2.0 / timestep * LVᵀmat(q2)' * Tmat() * Rmat(q3)' * Vᵀmat() * inertia * Vmat() * Lmat(q2)' * vector(q3)
        dynT = D2x + D1x
        dynR = D2q + D1q
        d = [dynT; dynR]

        # inputs
        barycenter_w = Dojo.vector_rotate(barycenter, q2) # TODO could be q3
        vc = v25 + Dojo.vector_rotate(Dojo.skew(-barycenter) * ϕ25, q2) # TODO could be q3
        vc_normal = vc' * impact_normal * impact_normal
        vc_tangential = vc - vc_normal

        F_impact = -opts.impact_damper * ψ * vc
        F_impact += opts.impact_spring * ψ * impact_normal
        F_friction = opts.sliding_friction * norm(F_impact) * coulomb_direction(vc_tangential)
        F_contact = F_impact + F_friction
        F2 = mass*gravity + F_contact

        τ_w = -Dojo.skew(F_contact) * barycenter_w
        τ2 = Dojo.vector_rotate(τ_w, inv(q2)) # TODO could be q3
        τ2 += -opts.rolling_drag * ψ * ϕ25
        τ2 += opts.rolling_friction * norm(F_impact) * coulomb_direction(ϕ25)

		d -= timestep * [F2; τ2]
        return d
    end

    vϕ = newton_solver(residual)
	v25 = vϕ[SUnitRange(1,3)]
	ϕ25 = vϕ[SUnitRange(4,6)]
	x3 = Dojo.next_position(x2, v25, timestep)
	q3 = Dojo.next_orientation(q2, ϕ25, timestep)
    z1 = [x3; v25; Dojo.vector(q3); ϕ25]
    return z1, ψ
end

function newton_solver(residual)
	x = zeros(6)
	for i = 1:10
		res = residual(x)
		(norm(res, Inf) < 1e-6) && break
		jac = FiniteDiff.finite_difference_jacobian(x -> residual(x), x)
		Δ = - jac \ res
		α = inegrator_line_search(x, res, Δ, residual)
		x += α * Δ
	end
	return x
end

function inegrator_line_search(x, res, Δ, residual)
    α = 1.0
    for i = 1:10
        (norm(residual(x + α * Δ), Inf) <= norm(res, Inf)) && break
        α *= 0.5
    end
    return α
end

function implicit_simulation(halfspace, soft, timestep, N, z0)
    z = [z0]
    ψ = []
    for i = 1:N
        @show i
        z0, ψ0 = implicit_integrator(halfspace, soft, timestep, z0)
        push!(z, z0)
        push!(ψ, ψ0)
    end
    return z, ψ
end



################################################################################
# Offline
################################################################################
halfspace_origin = [0,0,0.0]
normal = [0.0,0,1.0]
halfspace = HalfSpaceCollider(halfspace_origin, normal)
# soft = SoftCollider(nerf_object, mesh, N=5000)
soft = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]
@time collision(halfspace, soft)

set_light!(vis)
set_floor!(vis, color=RGBA(0.4,0.4,0.4,0.8), z=0.5, alt=0.0)
set_background!(vis)

timestep = 0.01
H = Int(floor(4.0/timestep))
mass = soft.mass
inertia = soft.inertia
gravity = [0,0,-9.81]

x2 = [0,0,1.0]
v15 = [0,2,3.0]
q2 = rand(4)
q2 ./= norm(q2)
q2 = Quaternion(q2...)
ϕ15 = [3.5,0,0.0]
z0 = [x2; v15; Dojo.vector(q2); ϕ15]

soft.options = ColliderOptions18(
    impact_damper=300.0,
    impact_spring=100.0,
    sliding_friction=0.2,
    coulomb_smoothing=1000.0,
    rolling_drag=0.5,
    rolling_friction=0.04,
)
@time Z0, Φ0 = implicit_simulation(halfspace, soft, timestep, H, z0)
plt = plot(layout=(3,1), xlabel="time")
plot!(plt[1,1], [i*h for i=0:H], [z[1] for z in Z0], linewidth=3.0, label="x", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[2] for z in Z0], linewidth=3.0, label="y", color=:blue)
plot!(plt[1,1], [i*h for i=0:H], [z[3] for z in Z0], linewidth=3.0, label="z", color=:blue)
plot!(plt[2,1], [i*h for i=0:H], [z[7] for z in Z0], linewidth=3.0, label="q1", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[8] for z in Z0], linewidth=3.0, label="q2", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[9] for z in Z0], linewidth=3.0, label="q3", color=:green)
plot!(plt[2,1], [i*h for i=0:H], [z[10] for z in Z0], linewidth=3.0, label="q4", color=:green)
plot!(plt[3,1], [i*h for i=1:H], Φ0, linewidth=3.0, label="ψ", color=:red)



build_mesh!(tight_mesh, vis, name=:tight_mesh, color=RGBA(0.6,0.8,0.9,1.0))
build_mesh!(mesh, vis, name=:mesh, color=RGBA(1,1,1,0.3))
# build_polytope!(hull, vis, name=:hull, color=RGBA(1,0.3,0.3,0.3))
animation = MeshCat.Animation(Int(floor(1/timestep)))
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



function set_floor!(vis::Visualizer;
	    x=20.0,
	    y=20.0,
	    z=0.1,
	    origin=[0,0,0.0],
		normal=[0,0,1.0],
	    color=RGBA(0.5,0.5,0.5,1.0),
	    tilepermeter=1.0,
	    imagename="tile.png",
	    axis::Bool=false,
	    grid::Bool=false)
    image = PngImage(joinpath(module_dir(), "assets", imagename))
    repeat = Int.(ceil.(tilepermeter * [x, y]))
    texture = Texture(image=image, wrap=(1,1), repeat=(repeat[1],repeat[2]))
    mat = MeshPhongMaterial(map=texture)
    (color != nothing) && (mat = MeshPhongMaterial(color=color))
    obj = HyperRectangle(Vec(-x/2, -y/2, -z), Vec(x, y, z))
    setobject!(vis[:floor], obj, mat)
	p = origin
	q = axes_pair_to_quaternion([0,0,1.], normal)
    settransform!(vis[:floor], MeshCat.compose(
		MeshCat.Translation(p...),
		MeshCat.LinearMap(q),
		))

    setvisible!(vis["/Axes"], axis)
    setvisible!(vis["/Grid"], grid)
    return nothing
end

set_floor!(vis, x=2.0, y=1.0, origin=[1,1,0.0], normal=[1,0,1.], grid=true)
MeshCat.MeshPhongMaterial()
x = 1.0
y = 1.0
z = 1.0
obj = HyperRectangle(Vec(0., 0, 0), Vec(x, y, z))
setobject!(vis[:floor], obj, mat)
settransform!(vis[:floor], MeshCat.Translation(-x/2, -y/2, -z))
