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

jr1.spring = 1e1 .* sones(3)# 1e4
jr1.damper = 1e1 .* sones(3)# 1e4
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
# Damper Jacobian
################################################################################

j0 = mech.eqconstraints[1]
jt0 = j0.constraints[1]
jr0 = j0.constraints[2]
origin0 = mech.origin
body0 = mech.bodies[2]
childid0 = 2
Δt0 = mech.Δt
damperforceb(jt0, origin0, body0, childid0, Δt0)
damperforceb(jr0, origin0, body0, childid0, Δt0)

damperforcea(joint::AbstractJoint, statea::State, stateb::State, Δt) = damperforcea(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., fullargssol(stateb)..., Δt)
damperforceb(joint::AbstractJoint, statea::State, stateb::State, Δt) = damperforceb(joint, posargsnext(statea, Δt)..., posargsnext(stateb, Δt)..., fullargssol(statea)..., fullargssol(stateb)..., Δt)
damperforceb(joint::AbstractJoint, stateb::State, Δt) = damperforceb(joint, posargsnext(stateb, Δt)..., fullargssol(stateb)..., Δt)

x2a0, q2a0 = posargsnext(statea, Δt)
x2b0, q2b0 = posargsnext(stateb, Δt)
x1a0, v1a0, q1a0, ω1a0 = fullargssol(statea)
x1b0, v1b0, q1b0, ω1b0 = fullargssol(stateb)

damperforcea(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)


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
