function getpendulum(; timestep::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
    spring = 0.0, damper = 0.0, spring_offset = szeros(1)) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    body1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_body1 = JointConstraint(Revolute(origin, body1,
        joint_axis; p2=p2, spring = spring, damper = damper, rotapply_springoffset = spring_offset,
        rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]))
    bodies = [body1]
    eqcs = [joint_between_origin_and_body1]

    mech = Mechanism(origin, bodies, eqcs, g = g, timestep = timestep, spring=spring, damper=damper)
    return mech
end

vis = Visualizer()
open(vis)
mech = getpendulum(timestep = 0.01, g = -9.81, spring = 0.0, damper = 0.0)
# mech.joints[1].λsol[2]
# reset!.(mech.joints)

ϕ1 = 0.4 * π
initialize!(mech, :pendulum, ϕ1 = ϕ1)
storage = simulate!(mech, 1.0, record = true, verbose = false)
visualize(mech, storage, vis=vis)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
Nb = length(collect(mech.bodies))
attjac = attitude_jacobian(data, Nb)

data = get_data(mech)
v15 = data[4:6]
sol = get_solution(mech)
v25 = sol[6:8]
norm(v15 - v25)

# IFT
set_entries!(mech)
datamat = full_data_matrix(mech, attjac = true)
datamat1 = full_data_matrix(mech, attjac = false)
datamat2 = full_data_matrix(mech, attjac = false) * attjac
plot(Gray.(attjac))

solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show cond(solmat)
@show rank(solmat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
fd_datamat1 = finitediff_data_matrix(mech, data, sol)
@test norm(fd_datamat + datamat, Inf) < 1e-7
@test norm(fd_datamat1 + datamat1, Inf) < 1e-7
@test norm(fd_datamat + datamat2, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))


fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(fd_solmat)))
plot(Gray.(abs.(solmat)))


fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
norm(fd_sensi - sensi, Inf)
norm(fd_sensi, Inf)

###############################################################################
# plot
###############################################################################

plot(hcat(Vector.(storage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')


sdf = get_sdf(mech, storage)
plot(hcat(sdf[1]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[2]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[3]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[4]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[5]...)')
plot(hcat(sdf[6]...)')
plot(hcat(sdf[7]...)')
plot(hcat(sdf[8]...)')
