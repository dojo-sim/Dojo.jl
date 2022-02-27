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
# using Plots
# using Random
# using MeshCat

# Open visualizer
vis=visualizer()
open(vis)

# Include new files
# include(joinpath(module_dir(), "examples", "loader.jl"))

mech = getmechanism(:box, timestep=0.01, g=-9.81, friction_coefficient=0.2, contact=false, mode=:box)#, contact_type = :nonlinear);
# mech = getmechanism(:box, timestep=0.01, g = -9.81, friction_coefficient = 0.2, contact = true, mode=:box, contact_type = :linear)
# mech = getmechanism(:box, timestep=0.01, g = -9.81, contact = true, mode=:box, contact_type = :impact)

Random.seed!(100)
ω = 10.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = 1.0 * [4, 1, 0.0]
initialize!(mech, :box, x = x, v = v, ω = ω)
storage = simulate!(mech, 1.3, record=true, verbose=false)

visualize(mech, storage, vis=vis)

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
# Set data
data = get_data(mech)
set_data!(mech, data)
sol = get_solution(mech)
Nb = length(mech.bodies)
attjac = attitude_jacobian(data, Nb)

# IFT
set_entries!(mech)
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show cond(solmat)
@show rank(solmat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

norm((datamat + fd_datamat)[1:6,:], Inf)
norm((datamat + fd_datamat)[1:6,1:6], Inf)
norm((datamat + fd_datamat)[1:6,7:12], Inf)
norm((datamat + fd_datamat)[1:6,13:18], Inf)
norm((datamat + fd_datamat)[7:12,:], Inf)
# norm((datamat + fd_datamat)[13:18,:], Inf)
norm((datamat + fd_datamat)[13:18,1:6], Inf)
# norm((datamat + fd_datamat)[13:18,7:12], Inf)
norm((datamat + fd_datamat)[13:18,13:18], Inf)

norm((datamat + fd_datamat)[1:6,7:12], Inf)
norm((datamat + fd_datamat)[1:3,7:12], Inf)
norm((datamat + fd_datamat)[4:6,7:12], Inf)
norm((datamat + fd_datamat)[4:6,7:9], Inf)
norm((datamat + fd_datamat)[4:6,10:12], Inf)

# norm((datamat + fd_datamat)[13:18,7:12], Inf)
# norm((datamat + fd_datamat)[13:15,7:12], Inf)
# norm((datamat + fd_datamat)[13:15,7:9], Inf)
norm((datamat + fd_datamat)[13:15,10:12], Inf)
norm((datamat + fd_datamat)[16:18,10:12], Inf)

# norm((datamat + fd_datamat)[4:6,7:9], Inf)
norm((datamat + fd_datamat)[13:15,7:9], Inf)

datamat[13:15,7:9] + fd_datamat[13:15,7:9]
datamat[13:15,7:9]
fd_datamat[13:15,7:9]


fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(fd_solmat)))
plot(Gray.(abs.(solmat)))

# plot(Gray.(abs.(1e13 * fd_solmat)))
# plot(Gray.(abs.(1e13 * solmat)))
# plot(Gray.(abs.(1e13 * (solmat + fd_solmat))))
#
# norm((solmat + fd_solmat)[1:3,:], Inf)
# norm((solmat + fd_solmat)[3:6,:], Inf)
# norm((solmat + fd_solmat)[4:6,1:3], Inf)
# norm((solmat + fd_solmat)[4:6,4:6], Inf)
# norm((solmat + fd_solmat)[4:6,7:8], Inf)
# norm((solmat + fd_solmat)[7:8,:], Inf)
#
# norm((solmat + fd_solmat)[4:6,4:6], Inf)
# fd_solmat[4:6,4:6]
# solmat[4:6,4:6]

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
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
