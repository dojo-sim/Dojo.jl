
vis = Visualizer()
open(vis)
render(vis)
include("env.jl")
include("initialize.jl")
# mech = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mech = get_rexhopper(timestep=0.01, gravity=0.0 * -9.81, model="rexhopper_damper",
    floating=true, contact=false, limits=false, spring=1.0, damper=0.5, contact_type=:nonlinear)
initialize!(mech, :rexhopper, x=[0,0,0.4])
# set_state!(mech, z0)
# z0 = get_maximal_state(mech)

loop_joint = get_joint_constraint(mech, :loop_joint)
constraint(mech, loop_joint)


function ctrl!(m,k)
    set_control!(m, 0sin(4k*m.timestep) * m.timestep * [szeros(6); sones(10)])
    return nothing
end
storage = simulate!(mech, 2.0, ctrl!, record=true, verbose=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)

mech.joints
mech.bodies

# convert_frames_to_video_and_gif("rexhopper_bounce_contact")
mech.bodies[1].state.x2
mech.bodies[1].state.q2

plot(hcat(get_sdf(mech, storage)[1][45:end]...)')
