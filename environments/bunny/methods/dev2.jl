using Dojo

vis = Visualizer()
open(vis)


mech = get_bunny(timestep=0.01)
mech.contacts[1].model.collision.collider.options = ColliderOptions()

initialize!(mech, :bunny,
    position=[0,0,0.6],
    # orientation=Quaternion(0,1,0,0.0,true),
    orientation=Quaternion(1,0,0,0.0,true),
    velocity=0*[0,0.5,5.0],
    angular_velocity=0*[0.5,10.0,3.0])

constraint(mech, mech.contacts[1])

@elapsed storage = simulate!(mech, 2.0,
    opts=SolverOptions(verbose=true, rtol=1e-4))
visualize(mech, storage, vis=vis)
mech.contacts[1]
constraint(mech, mech.contacts[1])


mech = get_bunny_sphere(timestep=0.01, gravity=-9.81)
mech.contacts[1].model.collision.collider.options = ColliderOptions()
mech.contacts[3].model.collision.collider.options = ColliderOptions()
initialize!(mech, :bunny_sphere,
    bunny_position=[0,0,0],
    bunny_velocity=[0,0,0],
    sphere_position=[0,4.0,0.4],
    sphere_velocity=[0,-5.0,0],
    )


# Main.@profiler
@elapsed storage = simulate!(mech, 5.0,
    opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)
mech.contacts[1]
mech.contacts[2]
mech.contacts[3]



mech = get_bunny_triumvirate(timestep=0.01, gravity=-9.81)
for i in [1,2,4,5]
    mech.contacts[i].model.collision.collider.options = ColliderOptions()
end

initialize!(mech, :bunny_triumvirate,
    positions=[[0,0,0.], [0,2,0.], [0,4,0.]],
    velocities=[[0,0,0.], [0,0,0.], [0,-5,0.]],
    )


# Main.@profiler
@elapsed storage = simulate!(mech, 5.0, opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)
mech.contacts[1]



collider = mech.contacts[1].model.collision.collider

particles_2 = zeros(Float32,10,3)
global const OSF_PATH = joinpath("/home/simon/research/repos/osf-pytorch")
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()
py"density_query"(collider.nerf_object, particles_2)
# py"density_query"(nerf_object, convert.(Float32, particles_2))









using FiniteDiff
using Plots

timestep = 0.01
model = mech.contacts[1].model
model.collider.options.impact_spring = 1e4
model.collider.options.impact_damper = 2e4
model.collider.options.sliding_drag = 3e0
model.collider.options.sliding_friction = 4e0
model.collider.options.rolling_drag = 1e0
model.collider.options.rolling_friction = 1e0
model.collider.options.coulomb_smoothing = 3e0
model.collider.options.coulomb_regularizer = 1e-0

vp = srand(3)
xp = srand(3)
qp = rand(4)
qp = Quaternion(qp ./ norm(qp)..., true)
ϕp = srand(3)

xc = srand(3)
vc = srand(3)
qc = rand(4)
qc = Quaternion(qc ./ norm(qc)..., true)
ϕc = srand(3)

x2p = next_position(xp, -vp, timestep)
q2p = next_orientation(qp, -ϕp, timestep)
x2c = next_position(xc, -vc, timestep)
q2c = next_orientation(qc, -ϕc, timestep)

function local_constraint(model, x2p, vp, q2p, ϕp, x2c, vc, q2c, ϕc, timestep)
    xp = next_position(x2p, vp, timestep)
    qp = next_orientation(q2p, ϕp, timestep)
    xc = next_position(x2c, vc, timestep)
    qc = next_orientation(q2c, ϕc, timestep)
    return constraint(model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
end

J0 = constraint_jacobian_velocity(:parent, model, xp, vp, qp, ϕp, xc, vc, qc, ϕc, timestep)
J1 = FiniteDiff.finite_difference_jacobian(
    vϕ -> local_constraint(model,
    x2p, vϕ[SUnitRange(1,3)], q2p, vϕ[SUnitRange(4,6)], x2c, vc, q2c, ϕc, timestep),
    [vp; ϕp])
norm(J0 - J1, Inf)
plot(Gray.(100abs.(Matrix(J0 - J1))))


v = srand(3)
J0 = 1/norm(v) * Diagonal(sones(3)) - v*v' ./ norm(v)^3
J1 = FiniteDiff.finite_difference_jacobian(v -> v/norm(v), v)
J0 - J1






















deps_folder = joinpath(module_dir(), "environments/bunny/deps")
inner_mesh_path = joinpath(deps_folder, "bunny_inner_mesh.obj")
outer_mesh_path = joinpath(deps_folder, "bunny_outer_mesh.obj")
collider = jldopen(joinpath(deps_folder, "soft_collider.jld2"))["soft"]
soft_body = SoftBody(collider, inner_mesh_path, outer_mesh_path)
[soft_body,]



mech.bodies[1].shape
Body(1.0, szeros(3,3), shape=mech.bodies[1].shape)
collider = mech.contacts[1].model.collider
bodies = [SoftBody(collider, name=:bunny, color=RGBA(1,1,1,0.1))]

# shape0 = bodies[1].shape
shape1 = Mesh(joinpath(deps_folder, "bunny_inner_mesh.obj"),
    position_offset = -collider.center_of_mass,
    color=RGBA(0.2,0.2,0.2,1.0))
shape2 = Mesh(joinpath(deps_folder, "bunny_outer_mesh.obj"),
    position_offset = -collider.center_of_mass,
    color=RGBA(0.9,0.9,0.9,0.3))
shape_vec = Vector{Shape{T}}([shape1, shape2])#, shape0]
shapes = Shapes(shape_vec)
bodies[1].shape = shapes



using Plots
using Dojo

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
light_blue = RGBA(0.4,0.4,1.0,0.8)

f(x) = x
ϕ(x) = x
L(x, ρ) = f(x) - ρ * log.(ϕ(x))
X = 0:0.0001:1
Y = -1:0.1:1
plt = plot(ylims=(-0.2,1.0), xlims=(-0.2,1.0), yticks=[0,1], xticks=[0,1], legend=:bottomright, size=(300,300))
plot!(plt, Y, f.(Y), linewidth=6.0, color=:black, label="f(x)")
plot!(rectangle(-0.2,2,0,-1), opacity=.6, color=:red, label="ϕ(x) < 0")
anim = @animate for i = 1:10
    plot!(plt, X, [1; L.(X[2:end], 10*0.5^i)], linewidth=5.0, color=light_blue, label=false)# label="ρ = 3e-1")
end
gif(anim, fps=10, "/home/simon/Downloads/anim.gif")

PHI = Vector(1:-0.005:0)
plt = plot(ylims=(-0.01,1.0), xlims=(-0.01,1.0), yticks=[0,1], xticks=[0,1], legend=:topright, size=(300,300), )
anim = @animate for i = 1:440
    # for (j,κ) ∈ enumerate([1e-1, 5e-2, 1e-2])
    for (j,κ) ∈ enumerate([5e-2,])
        PHI_i = PHI[1:min(i,length(PHI))]
        (i == 1) && plot!(plt, PHI_i, κ ./ PHI_i, color=light_blue, linewidth=1+1.5j, label="κ = $κ")
        plot!(plt, PHI_i, κ ./ PHI_i, color=light_blue, linewidth=1+1.5j, label=nothing)
    end
end

gif(anim, fps=40, "/home/simon/Downloads/force_solo.gif")

plot!(plt, PHI, 1e-2 ./ PHI, color=light_blue, linewidth=5.0, label="κ = 1e-2")
plot!(plt, PHI, 1e-3 ./ PHI, color=light_blue, linewidth=5.0, label="κ = 1e-3")
