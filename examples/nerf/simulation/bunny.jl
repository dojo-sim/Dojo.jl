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
# Export trajectories
################################################################################
using CSV
using DataFrames

px = [x[1] for x in storage.x[1]]
py = [x[2] for x in storage.x[1]]
pz = [x[3] for x in storage.x[1]]
qw = [vector(q)[1] for q in storage.q[1]]
qx = [vector(q)[2] for q in storage.q[1]]
qy = [vector(q)[3] for q in storage.q[1]]
qz = [vector(q)[4] for q in storage.q[1]]

plot(x1)
plot(x2)
plot(x3)

df = DataFrame(
    px = px,
    py = py,
    pz = pz,
    qw = qw,
    qx = qx,
    qy = qy,
    qz = qz,
    )
print(df)

CSV.write(joinpath(@__DIR__, "bunny_trajectory.csv"), df)
