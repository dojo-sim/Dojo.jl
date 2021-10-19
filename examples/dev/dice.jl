################################################################################
# Development
################################################################################
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

linmech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = false, conetype = :linear)
socmech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = false, conetype = :soc)

################################################################################
# Simulation
################################################################################

Random.seed!(100)
ω = (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = [1, 0.3, 0.2]

initialize!(linmech, :dice, x = x, v = v, ω = ω)
initialize!(socmech, :dice, x = x, v = v, ω = ω)
linstorage = simulate!(linmech, 0.5, record = true, solver = :mehrotra!)
socstorage = simulate!(socmech, 0.5, record = true, solver = :mehrotra!)


################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "dev", "diff_tools_control.jl"))
# Set data
Nb = length(socmech.bodies)
Random.seed!(10)
ndata = datadim(socmech, quat = true)
data = rand(ndata)*0.05
setdata!(socmech, data)
mehrotra!(socmech, opts = InteriorPointOptions(rtol = 1e-12, btol = 1e-12))
sol = getsolution(socmech)
attjac = attitudejacobian(data, Nb)

# IFT
datamat = full_data_matrix(socmech)
solmat = full_matrix(socmech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(socmech, data, sol) * attjac
@test norm(fd_datamat + datamat) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_solmat = finitediff_sol_matrix(socmech, data, sol)
@test norm(fd_solmat + solmat) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(socmech, data) * attjac
@test norm(fd_sensi - sensi) < 1e-8
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

###############################################################################
# plot
###############################################################################

plot(hcat(Vector.(linstorage.x[1])...)')
plot(hcat(Vector.(socstorage.x[1])...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in linstorage.q[1]]...)')
plot(hcat([[q.w, q.x, q.y, q.z] for q in socstorage.q[1]]...)')
plot(hcat(Vector.(linstorage.v[1])...)')
plot(hcat(Vector.(socstorage.v[1])...)')
plot(hcat(Vector.(linstorage.ω[1])...)')
plot(hcat(Vector.(socstorage.ω[1])...)')

sdf = get_sdf(linmech, linstorage)
sdf = get_sdf(socmech, socstorage)
plot(hcat(sdf[1]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[2]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[3]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[4]...)', ylims = (-0.01,0.01))
plot(hcat(sdf[5]...)')
plot(hcat(sdf[6]...)')
plot(hcat(sdf[7]...)')
plot(hcat(sdf[8]...)')
