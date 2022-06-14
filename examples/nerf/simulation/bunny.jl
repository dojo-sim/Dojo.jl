using Dojo
using Plots
using OSFLoader
using DojoEnvironments


vis = Visualizer()
render(vis)

################################################################################
# Simulate nerf
################################################################################
sliding_friction = 0.01
mech = get_mechanism(:nerf, nerf=:bunny, timestep=0.002,
    friction_coefficient=sliding_friction,
    collider_options=ColliderOptions(
        impact_damper=1e6,
        impact_spring=1e5,
        sliding_drag=0.01,
        sliding_friction=sliding_friction
        ))
mech.contacts[1].model.collision

initialize!(mech, :nerf,
    position=[0,0,0.6],
    orientation=Quaternion(normalize([1,0.6,0,0.5])...,true),
    # orientation=Quaternion(1,0,0,0.0,true),
    velocity=2*[1.0,1.0,0.5],
    angular_velocity=1*[0.5,1.0,1.0])

constraint(mech, mech.contacts[1])

@elapsed storage = simulate!(mech, 7.0,
    opts=SolverOptions(verbose=true, rtol=1e-4))
visualize(mech, storage, vis=vis)

render_static(vis)
open("/home/simon/Downloads/bunny_trajectory_$(sliding_friction).html", "w") do file
    write(file, static_html(vis))
end

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

# plot(px)
# plot(py)
# plot(pz)

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

CSV.write(joinpath(@__DIR__, "bunny_trajectory_$(sliding_friction).csv"), df)
