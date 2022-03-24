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

mech = get_bunny(timestep=0.01)
storage = simulate!(mech, 1.40, verbose=true)
visualize(mech, storage, vis=vis)

visualize(mech, storage, vis=vis)

normal = [0,0,1.0]
soft_constraint(mech.bodies[1], normal, contact_type=:soft)
SoftContact(mech.bodies[1], normal)
