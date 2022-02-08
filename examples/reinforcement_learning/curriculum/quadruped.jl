using Dojo
using IterativeLQR
using LinearAlgebra

gravity = -9.81
dt = 0.05
cf = 0.8
damper = 5.0
spring = 0.0
ρ0 = 1e-2
env = quadruped(
    mode=:min,
    dt=dt,
    gravity=gravity,
    cf=cf,
    damper=damper,
    spring=spring,
	infeasible_control=true,
	opts_step=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5),
    opts_grad=SolverOptions(rtol=ρ0, btol=ρ0, undercut=1.5)
	)
open(env.vis)


vis = Visualizer()
open(vis)
mech = get_quadruped(damper=0.3, spring=0spring, gravity=-9.81, body_contact=true,
	joint_limits=[100*[-0.5, -0.5, -0.5,],100*[ 0.5,  0.5,  0.5,]],
	limits=false)
initialize!(mech, :quadruped)
z0 = get_maximal_state(mech)
x0 = maximal_to_minimal(mech, z0)
x0[1:6]
# x0[7:8] .+= 0.5
x0[7] += 0.5
x0[9:10]
x0[11:12]
x0[13:14]
x0[15:16]
x0[17:18]

z1 = minimal_to_maximal(mech, x0)
x1 = maximal_to_minimal(mech, z1)
storage = generate_storage(mech, [z1])
visualize(mech, storage, vis=vis)



function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = 0.0*[szeros(6); sones(nu-6)] * mechanism.timestep
    set_control!(mechanism, u)
    return
end
storage = simulate!(mech, 4.0, ctrl!, record=true, verbose=false, opts=SolverOptions(btol=1e-5))
visualize(mech, storage, vis=vis, show_contact=true)




getfield.(mech.joints, :name)
# convert_frames_to_video_and_gif("quadruped_body_contact")
