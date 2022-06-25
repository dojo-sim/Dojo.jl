using Dojo
using Plots
using OSFLoader
using DojoEnvironments


vis = Visualizer()
render(vis)

################################################################################
# Simulate nerf
################################################################################
mech = get_mechanism(:nerf, nerf=:bunny, collider_options=ColliderOptions(), timestep=0.01)
mech.contacts[1].model.collision.options = ColliderOptions()

initialize!(mech, :nerf,
    position=[0,0,0.6],
    orientation=Quaternion(0,1,0,0.0,true),
    # orientation=Quaternion(1,0,0,0.0,true),
    velocity=1*[0,0.5,5.0],
    angular_velocity=1*[0.5,10.0,3.0])

constraint(mech, mech.contacts[1])

@elapsed storage = simulate!(mech, 7.0,
    opts=SolverOptions(verbose=true, rtol=1e-4))
visualize(mech, storage, vis=vis)

################################################################################
# Simulate nerf & sphere
################################################################################
mech = get_mechanism(:nerf_sphere, nerf=:bunny, collider_options=ColliderOptions(),
    mass=10.0,
    timestep=0.01,
    gravity=-9.81)
mech.contacts[1].model.collision.options = ColliderOptions()
mech.contacts[3].model.collision.options = ColliderOptions()

initialize!(mech, :nerf_sphere,
    nerf_position=[0,0,0],
    nerf_velocity=[0,0,0],
    sphere_position=[0.5,4.0,0.4],
    sphere_velocity=[0,-7.0,0],
    )
# Main.@profiler
@elapsed storage = simulate!(mech, 5.0,
    opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)

open(vis)
convert_frames_to_video_and_gif("sphere_bunny_simulation_corl")
