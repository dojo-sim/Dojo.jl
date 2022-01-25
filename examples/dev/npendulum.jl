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

mech = getmechanism(:npendulum, Δt = 0.01, g = -9.81, Nb = 2)
initialize!(mech, :npendulum, ϕ1 = 0.5)
initialize_simulation!(mech, true)
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]

storage = simulate!(mech, 1.11, record = true, solver = :mehrotra!)
visualize(mech, storage, vis = vis)

################################################################################
# Differentiation
################################################################################

# include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
Nb = length(mech.bodies)
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
attjac = attitude_jacobian(data, Nb)

# IFT
datamat = full_data_matrix(mech, attjac = true)
datamat0 = full_data_matrix(mech, attjac = true)
datamat1 = full_data_matrix(mech, attjac = false)
datamat2 = full_data_matrix(mech, attjac = false) * attjac

solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
fd_datamat0 = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
fd_datamat1 = finitediff_data_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_datamat0 + datamat0, Inf) < 1e-8
@test norm(fd_datamat1 + datamat1, Inf) < 1e-8
@test norm(fd_datamat0 + datamat2, Inf) < 1e-8
plot(Gray.(abs.(datamat1 + fd_datamat1)))
norm((datamat1 + fd_datamat1)[1:5,:], Inf)
norm((datamat1 + fd_datamat1)[6:10,:], Inf)
norm((datamat1 + fd_datamat1)[11:12,:], Inf)

norm((datamat1 + fd_datamat1)[6:10,1:13], Inf)
norm((datamat1 + fd_datamat1)[6:10,14:26], Inf)
norm((datamat1 + fd_datamat1)[6:10,27:28], Inf)

norm((datamat1 + fd_datamat1)[6:10,1:3], Inf)
norm((datamat1 + fd_datamat1)[6:10,4:6], Inf)
norm((datamat1 + fd_datamat1)[6:10,7:10], Inf)
norm((datamat1 + fd_datamat1)[6:10,11:13], Inf)

norm((datamat1 + fd_datamat1)[6:10,7:10], Inf)
(datamat1 + fd_datamat1)[6:10,7:10]
datamat1[6:10,    7:10]





-fd_datamat1[6:8,7:10]
(datamat1 + fd_datamat1)[6:8,7:10]
LVᵀmat(body1.state.q2[1])
(datamat1 + fd_datamat1)[6:8,7:10] * LVᵀmat(body1.state.q2[1])


plot(Gray.(abs.(fd_datamat1)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 1e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
