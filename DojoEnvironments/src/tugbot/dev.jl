
# ## Setup
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
Pkg.instantiate()
using Dojo 

# Examples
Pkg.activate(joinpath(@__DIR__, "..", "..", "examples")) 
Pkg.instantiate()
using StaticArrays

# include("methods/initialize.jl") 
# ##

vis = Visualizer()
open(vis)

mech = get_mechanism(:tugbot, gravity=-1.81, timestep=0.10)

function ctrl!(mechanism::Mechanism, k::Int; x_target=SVector{3}(0,3,3.0), kp=2.0, kd=2.0)
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

# Dojo.convert_frames_to_video_and_gif("tugbot")
