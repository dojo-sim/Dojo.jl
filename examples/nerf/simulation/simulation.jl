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
mech = get_mechanism(:nerf, nerf=:bluesoap, collider_options=ColliderOptions(), timestep=0.01)
mech = get_mechanism(:nerf, nerf=:halfsoap, collider_options=ColliderOptions(), timestep=0.01)
mech.contacts[1].model.collision.collider.options = ColliderOptions()

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
    timestep=0.01, gravity=-9.81)
mech = get_mechanism(:nerf_sphere, nerf=:bluesoap, collider_options=ColliderOptions(),
    timestep=0.01, gravity=-9.81)
mech = get_mechanism(:nerf_sphere, nerf=:halfsoap, collider_options=ColliderOptions(),
    timestep=0.01, gravity=-9.81)
mech.contacts[1].model.collision.collider.options = ColliderOptions()
mech.contacts[3].model.collision.collider.options = ColliderOptions()

initialize!(mech, :nerf_sphere,
    nerf_position=[0,0,0],
    nerf_velocity=[0,0,0],
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



################################################################################
# Simulate 2 nerfs & sphere
################################################################################
collider_options = ColliderOptions(
    impact_damper=3e5,
    impact_spring=3e4,
    sliding_friction=0.10,
    coulomb_smoothing=1e2)

mech = get_mechanism(:nerf_triumvirate, nerf=[:bunny, :bunny],
    collider_options=collider_options, timestep=0.01, gravity=-9.81)
# mech = get_mechanism(:nerf_triumvirate, nerf=[:bluesoap, :bluesoap],
    # collider_options=ColliderOptions(), timestep=0.01, gravity=-9.81)
# mech = get_mechanism(:nerf_triumvirate, nerf=[:halfsoap, :halfsoap],
    # collider_options=ColliderOptions(), timestep=0.01, gravity=-9.81)
for i in [1,2,4,5,6]
    mech.contacts[i].model.collision.options = ColliderOptions(
        impact_damper=3e5,
        impact_spring=3e4,
        sliding_friction=0.25)
end

contact = get_contact(mech, :contact_nerf1_nerf2)
contact.model.collision.options = ColliderOptions(
    impact_damper=1e5,
    impact_spring=3e4,
    sliding_drag=0.2,
    sliding_friction=0.1,
    torsional_drag=0.1,
    torsional_friction=0.01,
    rolling_drag=0.1,
    rolling_friction=0.01,
    coulomb_smoothing=10.0,
    coulomb_regularizer=1e-8)


initialize!(mech, :nerf_triumvirate,
    positions=[[0.1,0,-0.12], [-0.1,1.0,-0.12], [0,2.5,0.0]],
    velocities=[[0,0,0.], [0,0,0.], [0,-7,0.]],
    orientations=[one(Quaternion), one(Quaternion), one(Quaternion)]
    )

# Main.@profiler storage = simulate!(mech, 0.10, opts=SolverOptions(verbose=true, rtol=1e-4))
@elapsed storage = simulate!(mech, 1.5,
    opts=SolverOptions(verbose=true, rtol=3e-4, btol=3e-4))
visualize(mech, storage, vis=vis)

################################################################################
# Export
################################################################################
vis, anim = visualize(mech, storage, vis=vis, name=:white, color=RGBA(0.9,0.9,0.9,1.0))
vis, anim = visualize(mech, storage, vis=vis, animation=anim, name=:black, color=RGBA(0.2,0.2,0.2,1.0))
open(vis)
# convert_frames_to_video_and_gif("triumvirate_simulation")
