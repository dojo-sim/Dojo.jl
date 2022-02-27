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

linmech = getmechanism(:snake, Nb = 5, timestep=0.02, g = -9.81, friction_coefficient = 0.2, contact = false, contact_type = :linear_contact)
socmech = getmechanism(:snake, Nb = 5, timestep=0.02, g = -9.81, friction_coefficient = 0.2, contact = false, contact_type = :nonlinear)

x = [0,-1.0,0]
v = [1,.3,4]
ω = [.1,.8,0]
ϕ1 = π/2
initialize!(linmech, :snake, x = x, v = v, ω = ω, ϕ1 = ϕ1)
initialize!(socmech, :snake, x = x, v = v, ω = ω, ϕ1 = ϕ1)

@elapsed linstorage = simulate!(linmech, 1.5, record=true, solver = :mehrotra!)
@elapsed socstorage = simulate!(socmech, 1.5, record=true, solver = :mehrotra!)

visualize(linmech, linstorage, vis=vis)
visualize(socmech, socstorage, vis=vis)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
Nb = length(socmech.bodies)
data = get_data(socmech)
set_data!(socmech, data)
sol = get_solution(socmech)
attjac = attitude_jacobian(data, Nb)

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

fd_sensi = finitediff_sensitivity(socmech, data, δ = 1e-5, ϵb = 1e-14, ϵr = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 2e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))


q0 = UnitQuaternion(rand(4)...)
q1 = UnitQuaternion(rand(4)...)
p0 = rand(3)
r0 = vector_rotate(vector_rotate(p0, q1), q0)
r1 = vector_rotate(p0, q0*q1)

r0 - r1
