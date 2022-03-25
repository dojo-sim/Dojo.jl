using Dojo

vis = Visualizer()
open(vis)

################################################################################
# build NeRf
################################################################################
@pyinclude(joinpath(OSF_PATH, "extract_density_julia_cpu.py"))
nerf_object = py"generate_test_nerf"()
# SOFT = SoftCollider(nerf_object, N=5000)
# jldsave(joinpath(example_dir(), "results", "soft.jld2"), soft=SOFT)
SOFT = jldopen(joinpath(example_dir(), "results", "soft.jld2"))["soft"]

mech = get_bunny(SOFT, timestep=0.01)
mech.contacts[1].model.collider.options.impact_spring = 1e3
mech.contacts[1].model.collider.options.impact_damper = 3e3

VERBOSE = false
initialize!(mech, :bunny)
storage = simulate!(mech, 2.50, opts=SolverOptions(verbose=true, rtol=1e-4))
visualize(mech, storage, vis=vis)


impulse_map(mech, mech.contacts[1], mech.bodies[1])

mech.contacts[1].model.collider.options

normal = [0,0,1.0]
soft_constraint(mech.bodies[1], normal, contact_type=:soft)
SoftContact(mech.bodies[1], normal)
shape1 = Mesh(joinpath(module_dir(), "environments/bunny/deps/bunny_inner_mesh.obj"), color=RGBA(0.2,0.2,0.2,1.0))
