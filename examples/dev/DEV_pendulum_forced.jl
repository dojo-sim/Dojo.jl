# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))

# Open visualizer
# vis = Visualizer()
# open(vis)

# Build mechanism
mech = getmechanism(:pendulum, Δt = 0.01, g = -9.81)
initialize!(mech, :pendulum, ϕ1 = 0.7)

jointid = mech.eqconstraints[1].id
angles = zeros(1)
function controller!(mechanism, k)
    j1 = geteqconstraint(mechanism, jointid)
    θ1 = minimalCoordinates(mechanism, j1)[1]
    dθ1 = minimalVelocities(mechanism, j1)[1]
    u1 = (100.0*(angles[1]-θ1) + 5.0*(0-dθ1)) * mechanism.Δt
    setForce!(mechanism, j1, SA[u1])
    return
end

j1 = mech.eqconstraints[1]
jt1 = j1.constraints[1]
jr1 = j1.constraints[2]
j1.isdamper = true
j1.isspring = true

jr1.spring = 0.0 * sones(3)# 1e4
jr1.damper = 0.0 * sones(3)# 1e4
mech.eqconstraints[1].isdamper
mech.eqconstraints[1].constraints[2].damper

storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)
forcedstorage = simulate!(mech, 0.1, controller!, record = true, solver = :mehrotra!)
plot(hcat(Vector.(storage.x[1])...)')
plot(hcat(Vector.(forcedstorage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in forcedstorage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(forcedstorage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')
plot(hcat(Vector.(forcedstorage.ω[1])...)')

visualize(mech, storage, vis = vis)
# visualize(mech, forcedstorage, vis = vis)

################################################################################
# Differentiation
################################################################################

# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)

# IFT
datamat = full_data_matrix(deepcopy(mech))
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(deepcopy(mech), data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-8
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))
norm(fd_solmat + solmat, Inf)


fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 3e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
