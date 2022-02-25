
vis = Visualizer()
open(vis)
include("env.jl")
include("initialize.jl")
# mech = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mech = get_rexhopper(timestep=0.01, gravity= -0.99 * 9.81, model="rexhopper_no_wheel",
    floating=true, contact=true, limits=true, spring=0.0, damper=0.2, contact_type=:linear)

initialize!(mech, :rexhopper, x=[0,0,0.4])
z0 = get_maximal_state(mech)

# build_robot(mech, vis=vis, show_contact=true, color=RGBA(0.2, 0.2, 0.2, 1.0))
visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis, build=false)

function ctrl!(m,k)
    set_control!(m, 0*sin(4k*m.timestep) * m.timestep * [szeros(6); sones(10)])
    return nothing
end


VERBOSE = false
VERBOSE_MEHROTRA = false
# mech_failing = deepcopy(mech)
storage = simulate!(mech, 0.41, ctrl!, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))
visualize(mech, storage, vis=vis, show_contact=true, build=false)

mech.system.dimrow
mech.system.dimcol

mech.system.matrix_entries[1,8]
mech.bodies[1]
mech.system.matrix_entries[1,1]
mech.system.matrix_entries[2,2]
mech.system.matrix_entries[3,3]
mech.system.matrix_entries[4,4]
mech.system.matrix_entries[5,5]
mech.system.matrix_entries[6,6]
mech.system.matrix_entries[7,7]
mech.system.matrix_entries[7,8]
mech.system.matrix_entries[8,8]
mech.system.matrix_entries[9,9]
mech.system.matrix_entries[10,10]
mech.system.matrix_entries[11,11]
mech.system.matrix_entries[12,12]
mech.system.matrix_entries[13,13]
mech.system.matrix_entries[14,14]
mech.system.matrix_entries[15,15]

mech.joints[3].rotational.joint_limits
joint = mech.joints[3]

Diagonal(sones(3))

body = mech.bodies[1]
body.state.D[1:3,1:3] += rand(3,3)
# jldsave(joinpath(@__DIR__, "mech_failing.jld"), mech=mech_failing)
mech_failing = jldopen(joinpath(@__DIR__, "mech_failing.jld"))["mech"]
storage = simulate!(mech_failing, 0.01, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))

convert_frames_to_video_and_gif("rexhopper_new_failure")

mech.bodies[1].state.D



function control!(mechanism, k; u=0.1)
    for (i, joint) in enumerate(mechanism.joints)
        nu = Dojo.control_dimension(joint, ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        Dojo.set_input!(joint, su)
    end
    return
end

# mechanism
mech = get_mechanism(:pendulum)
initialize!(mech, :pendulum)

mech.bodies[1].state.ϕsol[2] = 0.5*srand(3)
mech.bodies[1].impulses[2] = 0.5*sones(1)
mech.bodies[1].impulses_dual[2] = 0.1*sones(1)
# simulate
# storage = simulate!(mech, 0.1, control!,
    # record=true, opts=SolverOptions(rtol=1e-4, btol=1e-4, verbose=true))

# Set data
Nb = length(mech.bodies)
data = Dojo.get_data0(mech)
Dojo.set_data0!(mech, data)
sol = Dojo.get_solution0(mech)
attjac = Dojo.attitude_jacobian(data, Nb)

# IFT
set_entries!(mech)
solmat = Dojo.full_matrix(mech.system)
# finite diff
fd_solmat = finitediff_sol_matrix(mech, data, sol, δ=1.0e-5, verbose=true)
norm(fd_solmat + solmat, Inf) < 1e-4

solmat
fd_solmat
solmat + fd_solmat
plot(Gray.(1e8*abs.(solmat+fd_solmat)))

plot(Gray.(solmat))
plot(Gray.(fd_solmat))


plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))
