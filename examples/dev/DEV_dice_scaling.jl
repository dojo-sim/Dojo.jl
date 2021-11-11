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
include(joinpath(module_dir(), "examples", "loader.jl"))

<<<<<<< HEAD
mech = getmechanism(:dice, Δt = 1.0, g = -1.0, cf = 1.0, contact = false, conetype = :contact)
=======
mech = getmechanism(:dice, Δt = 0.01, g = -9.81, cf = 0.2, contact = true, conetype = :soc)
>>>>>>> 5675fe2e79622be7440f9fd2c6e43ea0f2940ce1

Random.seed!(100)
ω = 0.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = 0.0 * [1, 0.3, 0.2]
initialize!(mech, :dice, x = x, v = v, ω = ω)
<<<<<<< HEAD
storage = simulate!(mech, 0.1, record = false, solver = :mehrotra!)
@show data
=======
storage = simulate!(mech, 0.01, record = true, solver = :mehrotra!)

>>>>>>> 5675fe2e79622be7440f9fd2c6e43ea0f2940ce1
################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools_control_contact.jl"))
# Set data
Nb = length(mech.bodies)
# Random.seed!(10)
# ndata = datadim(mech, quat = true)
# data = rand(ndata)*0.05
data = getdata(mech)
setdata!(mech, data)
mehrotra!(mech, opts = InteriorPointOptions(rtol = 1e-6, btol = 1e-1, undercut=1.2, verbose=true))
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

# IFT
setentries!(mech)
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show cond(solmat)
@show rank(solmat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
<<<<<<< HEAD
@test norm(fd_datamat + datamat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))
=======
@test norm(fd_datamat + datamat, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))
>>>>>>> 5675fe2e79622be7440f9fd2c6e43ea0f2940ce1

fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

<<<<<<< HEAD
fd_sensi_ift = -fd_solmat \ fd_datamat

@show data
fd_sensi = finitediff_sensitivity(mech, data, δ=1e-5, ϵr=1.0e-12, ϵb=1.0e-1) * attjac
@test norm(fd_sensi - sensi, Inf) / norm(fd_sensi, Inf)  < 1e-3
@test norm(fd_sensi - fd_sensi_ift, Inf)  < 1e-3
@test norm(sensi - fd_sensi_ift, Inf)  < 1e-3
norm(fd_sensi - sensi, Inf)
=======
fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
>>>>>>> 5675fe2e79622be7440f9fd2c6e43ea0f2940ce1
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
