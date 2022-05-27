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
using RobotVisualizer
using Quaternions

# include("methods/initialize.jl")
# ##

# ## visualizer
vis = Visualizer()
open(vis)

mech = DojoEnvironments.get_mechanism(:tugbot; gravity=-1.81, timestep=0.10, radius=0.32)

function ctrl!(mechanism::Mechanism, k::Int; x_target=SVector{3}(0,3,3.0), kp=1.5, kd=2.0)
    dt = mechanism.timestep
    drone_body = get_body(mechanism, :drone)
    drone_joint = get_joint(mechanism, :drone_joint)
    x = Dojo.current_position(drone_body.state)
    v, ω = Dojo.current_velocity(drone_body.state)
    θ = mrp(current_orientation(mech.bodies[1].state))
    u_gravity = -drone_body.mass * mechanism.gravity * dt
    u_tra = u_gravity + kp*(x_target - x)*dt - kd * dt .* v
    u_rot = -kp/5 * θ - kd/10 * dt .* ω
    set_input!(drone_joint, [u_tra; u_rot])
    return nothing
end

initialize!(mech, :tugbot, drone_position=[1,0,0])

storage = simulate!(mech, 5.2, ctrl!, record=true, verbose=true)
vis, anim = visualize(mech, storage, vis=vis, show_contact=false)

# add drone mesh
tugbot_visualizer!(vis, mech)

# find the trajectories of the two attach points
contact = get_contact(mech, :body_body)
parent_traj = []
child_traj = []
for i = 1:52
    z = get_maximal_state(storage, i)
    set_maximal_state!(mech, z)
    x_parent, x_child = string_location(mech, contact.model.collision)
    push!(parent_traj, Vector(x_parent))
    push!(child_traj, Vector(x_child))
end

# add rope to visualizer
N = 50
build_rope(vis, N=N, name=:rope2, rope_radius=0.01)
set_loose_rope(vis, [1,-1,1.0], [1,1,1.0], rope_length=10, N=N, min_altitude=-0.0, name=:rope2)

vis, anim = animate_loose_rope(vis, parent_traj, child_traj, rope_length=2.0, N=N,
    anim=anim, name=:rope2, min_altitude=0.02)

convert_frames_to_video_and_gif("ballincup")
