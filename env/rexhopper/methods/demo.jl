
vis = Visualizer()
open(vis)
include("env.jl")
include("initialize.jl")
# mech = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mech = get_rexhopper(timestep=0.02, gravity= -1.0 * 9.81, model="rexhopper_no_wheel0",
    floating=true, contact=true, limits=true, spring=0.0, damper=0.5, contact_type=:linear)

initialize!(mech, :rexhopper, x=[0,0,0.4])
# set_state!(mech, z0)
z0 = get_maximal_state(mech)
visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis)

function ctrl!(m,k)
    set_input!(m, 4.0 * sin(4k*m.timestep) * m.timestep * [szeros(6); sones(10)])
    return nothing
end
storage = simulate!(mech, 10.0, ctrl!, record=true, verbose=false,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=false))
visualize(mech, storage, vis=vis, show_contact=true)


# convert_frames_to_video_and_gif("rexhopper_bounce_contact")
mech.bodies[1].state.x2
mech.bodies[1].state.q2

plot(hcat(get_sdf(mech, storage)[1][45:end]...)')

build_robot(mech, vis=vis, show_contact=true, name=:hopper, color=RGBA(0.4, 0.4, 0.6, 1.0))
@elapsed visualize(mech, storage, vis=vis, show_contact=true, build=true, name=:hopper)
@elapsed visualize(mech, storage, vis=vis, show_contact=true, build=false, name=:hopper)


## quadruped 
mech = get_quadruped(timestep=0.01, gravity= -9.81, contact=true, body_contact=true,
    limits=true, spring=0.0, damper=0.5)
initialize!(mech, :quadruped, tran=[0,0,0.2], rot=[0.003,0.0003,0.])
z0 = get_maximal_state(mech)
visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis, build=true);
storage = simulate!(mech, 0.41, record=true, abort_upon_failure=false,
    opts=SolverOptions(rtol=1e-5, btol=1e-5, undercut=5.0, no_progress_undercut=10.0, verbose=true));
visualize(mech, storage, vis=vis, show_contact=true, build=true)

cond(full_matrix(mech.system))