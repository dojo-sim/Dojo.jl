
using Dojo
global REG = 1.0e-6


vis= Visualizer()
open(vis)
include("env.jl")
include("initialize.jl")
# mechanism = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mechanism = get_rexhopper(timestep=0.01, gravity= -0.99 * 9.81, model="rexhopper_no_wheel",
    floating=true, contact_foot=true, limits=true, springs=0, damper=0.5, contact_type=:linear)

q0 = [1,0.5,0,0]
q0 = Quaternion(q0 ./ norm(q0)...)
initialize!(mechanism, :rexhopper, body_position=[0,0,0.4], body_orientation=[0.5,0.8,0])
z0 = get_maximal_state(mechanism)

storage = simulate!(mechanism, 0.40, record=true, abort_upon_failure=false,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5, verbose=true))
visualize(mechanism, storage, vis=vis, show_contact=true, build=true)



# build_robot(mechanism, vis=vis, show_contact=true, color=RGBA(0.2, 0.2, 0.2, 1))
visualize(mechanism, generate_storage(mechanism, [z0]), show_contact=true, vis=vis, build=false)

function ctrl!(m,k)
    set_input!(m, 0*sin(4k) * [szeros(6); sones(10)])
    return nothing
end


verbose=false
VERBOSE_MEHROTRA = false
# mech_failing = deepcopy(mechanism)
storage = simulate!(mechanism, 0.41, ctrl!, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5, verbose=true))
visualize(mechanism, storage, vis=vis, show_contact=true, build=false)

mechanism.system.dimrow
mechanism.system.dimcol

mechanism.system.matrix_entries[1,8]
mechanism.bodies[1]
mechanism.system.matrix_entries[1,1]
mechanism.system.matrix_entries[2,2]
mechanism.system.matrix_entries[3,3]
mechanism.system.matrix_entries[4,4]
mechanism.system.matrix_entries[5,5]
mechanism.system.matrix_entries[6,6]
mechanism.system.matrix_entries[7,7]
mechanism.system.matrix_entries[7,8]
mechanism.system.matrix_entries[8,8]
mechanism.system.matrix_entries[9,9]
mechanism.system.matrix_entries[10,10]
mechanism.system.matrix_entries[11,11]
mechanism.system.matrix_entries[12,12]
mechanism.system.matrix_entries[13,13]
mechanism.system.matrix_entries[14,14]
mechanism.system.matrix_entries[15,15]

mechanism.joints[3].rotational.joint_limits
joint = mechanism.joints[3]

Diagonal(sones(3))

body = mechanism.bodies[1]
body.state.D[1:3,1:3] += rand(3,3)
# jldsave(joinpath(@__DIR__, "mech_failing.jld"), mech=mech_failing)
mech_failing = jldopen(joinpath(@__DIR__, "mech_failing.jld"))["mech"]
storage = simulate!(mech_failing, 0.01, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5, verbose=true))

convert_frames_to_video_and_gif("rexhopper_new_failure")

mechanism.bodies[1].state.D



function control!(mechanism, k; u=0.1)
    for (i, joint) in enumerate(mechanism.joints)
        nu = Dojo.control_dimension(joint, ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        Dojo.set_input!(joint, su)
    end
    return
end

# mechanism
mechanism = get_mechanism(:pendulum)
initialize!(mechanism, :pendulum)

mechanism.bodies[1].state.ωsol[2] = 0.5*srand(3)
mechanism.bodies[1].impulses[2] = 0.5*sones(1)
mechanism.bodies[1].impulses_dual[2] = 0.1*sones(1)
# simulate
# storage = simulate!(mechanism, 0.1, control!,
    # record=true, opts=SolverOptions(rtol=1e-4, btol=1e-4, verbose=true))

# Set data
Nb = length(mechanism.bodies)
data = Dojo.get_data0(mechanism)
Dojo.set_data0!(mechanism, data)
sol = Dojo.get_solution0(mechanism)
attjac = Dojo.attitude_jacobian(data, Nb)

# IFT
set_entries!(mechanism)
solmat = Dojo.full_matrix(mechanism.system)
# finite diff
fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ=1.0e-5, verbose=true)
norm(fd_solmat + solmat, Inf) < 1e-4

solmat
fd_solmat
solmat + fd_solmat
plot(Gray.(1e8*abs.(solmat+fd_solmat)))

plot(Gray.(solmat))
plot(Gray.(fd_solmat))


plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))




mechanism = get_quadruped(timestep=0.01, gravity= -9.81, contact_feet=true, contact_body=true,
    limits=true, springs=0, damper=0.2)

initialize!(mechanism, :quadruped, tran=[0,0,0.2], rot=[0.003,0.0003,0.])
z0 = get_maximal_state(mechanism)

visualize(mechanism, generate_storage(mechanism, [z0]), show_contact=true, vis=vis, build=true)

storage = simulate!(mechanism, 0.40, record=true, abort_upon_failure=false,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5, verbose=true))
visualize(mechanism, storage, vis=vis, show_contact=true, build=true)

sol0 = get_solution(mechanism)
data0 = get_data(mechanism)

scn(residual_violation(mechanism), digits=0)
residual_violation(mechanism)

set_entries!(mechanism)
fv0 = full_vector(mechanism.system)
fM0 = full_matrix(mechanism.system)
Δv0 = fM0 \ fv0

ldu_factorization!(mechanism.system)    # factorize system, modifies the matrix in place
ldu_backsubstitution!(mechanism.system) # solve system, modifies the vector in place
Δv1 = full_vector(mechanism.system)



mechanism.system.vector_entries
res = norm(fv0, Inf)
Δ0 = norm(Δv0, Inf)

using Plots
plot(Gray.(fM0))

findmax(fv0)
fM0[:,111]
plot(Gray.(1e1abs.(reshape(Δv0, (6,47)))))
plot(Gray.(1e1abs.(reshape(Δv1, (6,47)))))
plot(Gray.(1e10abs.(reshape(Δv0 - Δv1, (6,47)))))
norm(Δv0 - Δv1, Inf)

plot(Gray.(1e1abs.(reshape(fv0, (6,47)))))
plot(Gray.(1e1abs.(reshape(fM0[:,111], (6,47)))))

fv0[111]
fM0[111,111]
Δv0[111]

dc0 = mechanism.system.dimcol
[sum(dc0[1:i-1])+1:sum(dc0[1:i]) for i in 1:length(dc0)][14]

getfield.(mechanism.bodies, :id)
