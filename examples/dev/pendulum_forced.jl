# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

# Load packages
using Plots
using Random
using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

mech = getmechanism(:pendulum, Δt = 0.01, g = -9.81)
initialize!(mech, :pendulum, ϕ1 = 0.7)

jointid = mech.eqconstraints[1].id
angles = zeros(1)
function controller!(mechanism, k)
    j1 = get_joint_constraint(mechanism, jointid)
    θ1 = minimal_coordinates(mechanism, j1)[1]
    dθ1 = minimal_velocities(mechanism, j1)[1]
    u1 = (100.0*(angles[1]-θ1) + 5.0*(0-dθ1)) * mechanism.Δt
    set_input!(j1, SA[u1])
    return
end

storage = simulate!(mech, 4.0, record = true, solver = :mehrotra!)
forcedstorage = simulate!(mech, 4.0, controller!, record = true, solver = :mehrotra!)
plot(hcat(Vector.(storage.x[1])...)')
plot(hcat(Vector.(forcedstorage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in forcedstorage.q[1]]...)')
plot(hcat(Vector.(storage.v[1])...)')
plot(hcat(Vector.(forcedstorage.v[1])...)')
plot(hcat(Vector.(storage.ω[1])...)')
plot(hcat(Vector.(forcedstorage.ω[1])...)')

visualize(mech, storage, vis = vis)
visualize(mech, forcedstorage, vis = vis)

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

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 3e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
