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
vis=visualizer()
open(vis)

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

mech = getmechanism(:atlas, timestep=0.01, g = -9.81, friction_coefficient = 0.5, damper = 50.0, spring=1.0, contact = true)

initialize!(mech, :atlas, tran = [0,0,0.5], rot = [0.01,0.05, 0.0])
storage = simulate!(mech, 1.0, record=true, opts = SolverOptions(btol = 1e-6))
visualize(mech, storage, vis=vis)

# Set data
Nb = length(mech.bodies)
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
attjac = attitude_jacobian(data, Nb)

# IFT
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-6
plot(Gray.(abs.(1e10 .* datamat)))
plot(Gray.(abs.(1e10 .* fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(1e10 * fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 8e-3
plot(Gray.(1e10 .* sensi))
plot(Gray.(fd_sensi))


test_solmat(:atlas, ϵ = 1e-8)
test_datamat(:atlas, ϵ = 1e-6)
test_sensitivity(:atlas, ϵ = 8e-3)

test_solmat(:quadruped, ϵ = 1e-8)
test_datamat(:quadruped, ϵ = 1e-6)
test_sensitivity(:quadruped, ϵ = 8e-3)
