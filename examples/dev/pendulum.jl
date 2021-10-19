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
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))

# mech = getmechanism(:pendulum, Δt = 0.01, g = -9.81)
mech = getmechanism(:npendulum, Δt = 0.01, g = -9.81, Nlink = 1)
# initialize!(mech, :pendulum, ϕ1 = 0.7)
initialize!(mech, :npendulum, ϕ1 = 0.0)
storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "dev", "diff_tools_control.jl"))
# Set data
Nb = length(mech.bodies)
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

resvec = full_vector(mech.system)
norm(resvec, Inf)
# IFT
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
sensi2 = sensitivities(mech, sol, data)

@test norm(sensi - sensi2, Inf) < 1.0e-8
v0 = rand(13)
@test norm(jvp(mech, sol, data, v0) - sensi * v0, Inf) < 1.0e-8

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac

@test norm(fd_datamat + datamat, Inf) < 1e-8
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵ=1.0e-12) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 1e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

norm(fd_sensi - sensi, Inf) / norm(fd_sensi, Inf)
norm(fd_sensi - sensi, Inf)
