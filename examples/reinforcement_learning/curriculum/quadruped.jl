using Dojo
using IterativeLQR
using LinearAlgebra

vis = Visualizer()
open(vis)

# z0 = get_maximal_state(mech)
# x0 = maximal_to_minimal(mech, z0)
# x0[1:6]
# x0[7:12]
# # x0[13:14] .= -0.50 # RLhip -0.50, +0.50
# # x0[15:16] .= 1.5 # RLthigh -0.50, +1.50
# # x0[17:18] .= -2.50 # RLcalf -2.50, -1.00
# z1 = minimal_to_maximal(mech, x0)
# x1 = maximal_to_minimal(mech, z1)
# storage = generate_storage(mech, [z1])
# visualize(mech, storage, vis=vis)



gravity = -0*9.81
timestep = 0.01
friction_coefficient = 0.8
damper = 0.0
spring = 0.0
ρ0 = 1e-2

mech = get_quadruped(timestep=timestep, damper=damper, spring=spring,
	gravity=gravity, body_contact=true,
	joint_limits=[[-0.5, -0.5, -2.5,],[ 0.5,  1.5, -1.0,]],
	limits=true)
initialize!(mech, :quadruped)
function ctrl!(mechanism, k)
    nu = control_dimension(mechanism)
    u = -1*[szeros(6); sones(nu-6)] * mechanism.timestep
    set_control!(mechanism, u)
    return
end
storage = simulate!(mech, 2.0, ctrl!, record=true, verbose=true, opts=SolverOptions(btol=1e-5, verbose=true))
visualize(mech, storage, vis=vis, show_contact=true)
