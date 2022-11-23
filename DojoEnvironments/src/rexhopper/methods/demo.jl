
using Dojo
global REG = 1.0e-6


vis= Visualizer()
open(vis)
include("env.jl")
include("initialize.jl")
# mech = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mech = get_rexhopper(timestep=0.01, gravity= -0.99 * 9.81, model="rexhopper_no_wheel",
    floating=true, contact_foot=true, limits=true, spring=0.0, damper=0.5, contact_type=:linear)

q0 = [1,0.5,0,0]
q0 = Quaternion(q0 ./ norm(q0)...)
initialize!(mech, :rexhopper, body_position=[0,0,0.4], body_orientation=[0.5,0.8,0.0])
z0 = get_maximal_state(mech)

storage = simulate!(mech, 0.40, record=true, abort_upon_failure=false,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))
visualize(mech, storage, vis=vis, show_contact=true, build=true)



# build_robot(mech, vis=vis, show_contact=true, color=RGBA(0.2, 0.2, 0.2, 1.0))
visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis, build=false)

function ctrl!(m,k)
    set_input!(m, 0*sin(4k) * [szeros(6); sones(10)])
    return nothing
end


verbose=false
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

mech.bodies[1].state.ωsol[2] = 0.5*srand(3)
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




mech = get_quadruped(timestep=0.01, gravity= -9.81, contact_feet=true, contact_body=true,
    limits=true, spring=0.0, damper=0.2)

initialize!(mech, :quadruped, tran=[0,0,0.2], rot=[0.003,0.0003,0.])
z0 = get_maximal_state(mech)

visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis, build=true)

storage = simulate!(mech, 0.40, record=true, abort_upon_failure=false,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))
visualize(mech, storage, vis=vis, show_contact=true, build=true)

sol0 = get_solution(mech)
data0 = get_data(mech)

scn(residual_violation(mech), digits=0)
residual_violation(mech)

set_entries!(mech)
fv0 = full_vector(mech.system)
fM0 = full_matrix(mech.system)
Δv0 = fM0 \ fv0

ldu_factorization!(mech.system)    # factorize system, modifies the matrix in place
ldu_backsubstitution!(mech.system) # solve system, modifies the vector in place
Δv1 = full_vector(mech.system)



mech.system.vector_entries
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

dc0 = mech.system.dimcol
[sum(dc0[1:i-1])+1:sum(dc0[1:i]) for i in 1:length(dc0)][14]

getfield.(mech.bodies, :id)
