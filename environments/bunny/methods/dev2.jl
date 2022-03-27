using Dojo

vis = Visualizer()
open(vis)



mech = get_bunny(timestep=0.01)
mech.contacts[1].model.collision.collider.options.impact_spring = 3e4
mech.contacts[1].model.collision.collider.options.impact_damper = 3e5
mech.contacts[1].model.collision.collider.options.sliding_drag = 0.10
mech.contacts[1].model.collision.collider.options.sliding_friction = 0.2
mech.contacts[1].model.collision.collider.options.rolling_drag = 0.0
mech.contacts[1].model.collision.collider.options.rolling_friction = 0.01
mech.contacts[1].model.collision.collider.options.coulomb_smoothing = 3e+1
mech.contacts[1].model.collision.collider.options.coulomb_regularizer = 1e-3

mech.contacts[1]
mech.bodies[1].inertia = 0.1*Diagonal(sones(3))

# mech.origin.state.x1 = SVector{3}(0,0,0.0)
# mech.origin.state.x2 = SVector{3}(0,0,0.0)
# q0 = [1,1,0,0.]
# q0 ./= norm(q0)
# mech.origin.state.q1 = Quaternion(q0...,true)
# mech.origin.state.q2 = Quaternion(q0...,true)

VERBOSE = false
initialize!(mech, :bunny,
    position=[0,0,0.3],
    orientation=Quaternion(0,1,0,0.0,true),
    velocity=[0,0.5,5.0],
    angular_velocity=[0.5,10.0,3.0])

constraint(mech, mech.contacts[1])

storage = simulate!(mech, 6.0, opts=SolverOptions(verbose=true, rtol=1e-4, max_iter=50))
visualize(mech, storage, vis=vis)



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
