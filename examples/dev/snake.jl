using ConstrainedDynamics
using ConstrainedDynamicsVis

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

# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))

linmech = getmechanism(:snake, Nlink = 5, Δt = 0.01, g = -9.81, cf = 0.2, contact = false, conetype = :linear)
socmech = getmechanism(:snake, Nlink = 5, Δt = 0.01, g = -9.81, cf = 0.2, contact = false, conetype = :soc)

x = [0,-0.5,0]
v = [1,.3,4]
ω = [.1,.8,0]
ϕ1 = π/2
initialize!(linmech, :snake, x = x, v = v, ω = ω, ϕ1 = ϕ1)
initialize!(socmech, :snake, x = x, v = v, ω = ω, ϕ1 = ϕ1)

linstorage = simulate!(linmech, 0.5, record = true, solver = :mehrotra!)
socstorage = simulate!(socmech, 0.5, record = true, solver = :mehrotra!)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "dev", "diff_tools_control.jl"))
# Set data
Nb = length(socmech.bodies)
data = getdata(socmech)
setdata!(socmech, data)
sol = getsolution(socmech)
attjac = attitudejacobian(data, Nb)

# IFT
datamat = full_data_matrix(socmech)
solmat = full_matrix(socmech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(socmech, data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-8
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(socmech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(socmech, data, δ = 1e-5, ϵ = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 2e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
