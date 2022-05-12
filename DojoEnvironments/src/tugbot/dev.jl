# ## Setup
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
Pkg.instantiate()
using Dojo

# Examples
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "examples"))
Pkg.instantiate()
using StaticArrays
using DojoEnvironments

# include("methods/initialize.jl")
# ##

# ## visualizer
vis = Visualizer()
open(vis)

mech = DojoEnvironments.get_mechanism(:tugbot, gravity=-1.81, timestep=0.10, radius=0.30)
mech.bodies[1].inertia *= 10.0

function ctrl!(mechanism::Mechanism, k::Int; x_target=SVector{3}(0,3,3.0), kp=1.5, kd=2.0)
    dt = mechanism.timestep
    drone_body = get_body(mechanism, :drone)
    drone_joint = get_joint(mechanism, :drone_joint)
    x = Dojo.current_position(drone_body.state)
    v = Dojo.current_velocity(drone_body.state)[1]
    u_gravity = -drone_body.mass * mechanism.gravity * dt
    u_tra = u_gravity + kp*(x_target - x)*dt - kd * dt .* v
    u_rot = szeros(3)
    set_input!(drone_joint, [u_tra; u_rot])
    return nothing
end

initialize!(mech, :tugbot, drone_position=[1,0,0])

storage = simulate!(mech, 5.2, ctrl!, record=true, verbose=true)
visualize(mech, storage, vis=vis, show_contact=false)

function tugbot_visualizer!(vis::Visualizer, mech::Mechanism;)
    # image = PngImage(image1)
    # texture = Texture(image=image);
    material = MeshPhongMaterial(color=Colors.RGBA(0.0, 0.0, 0.0, 1.0))
    drone = MeshFileObject(joinpath(dirname(pathof(DojoEnvironments)), "tugbot/deps/drone.obj"))
    setobject!(vis["robot"]["bodies"]["drone__id_1"][:mesh], drone)
    settransform!(vis["robot"]["bodies"]["drone__id_1"][:mesh], MeshCat.LinearMap(RotX(π/2)))
    # image = PngImage(image2)
    # texture = Texture(image=image);
    # material = MeshLambertMaterial(map=texture);

    # setobject!(vis["bodies"]["body:2"], convert_shape(mech.bodies[4].shape), material)
end


tugbot_visualizer!(vis, mech)
# Dojo.convert_frames_to_video_and_gif("tugbot")
