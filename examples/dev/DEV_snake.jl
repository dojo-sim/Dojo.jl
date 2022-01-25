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


mech = getmechanism(:snake, Nb=5, Δt=0.01, g=-9.81, cf=0.0, contact=true, contact_type=:contact)

x = [0,-0.5,.1]
v = 0.1*[1,.3,4]
ω = 1.5*[.1,.8,0]
ϕ1 = π/1.5
initialize!(mech, :snake, x = x, v = v, ω = ω, ϕ1 = ϕ1)

@elapsed storage = simulate!(mech, .28, record = true, solver = :mehrotra!)

visualize(mech, storage, vis = vis)

################################################################################
# Differentiation
################################################################################

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
norm(fd_datamat + datamat, Inf)

# norm((fd_datamat + datamat)[1:10,1:44], Inf)
# norm((fd_datamat + datamat)[11:16,1:44], Inf)
# norm((fd_datamat + datamat)[17:22,1:44], Inf)
# norm((fd_datamat + datamat)[23:28,1:44], Inf)
# norm((fd_datamat + datamat)[29:end,1:44], Inf)
#
# norm((fd_datamat + datamat)[17:22,13:24], Inf)
# norm((fd_datamat + datamat)[23:28,25:36], Inf)
# norm((fd_datamat + datamat)[23:28,25:36], Inf)



norm((fd_datamat + datamat)[1:5,1:31], Inf)
norm((fd_datamat + datamat)[6:17,1:31], Inf)

norm((fd_datamat + datamat)[6:17,1:12], Inf)
norm((fd_datamat + datamat)[6:17,13:24], Inf)
norm((fd_datamat + datamat)[6:17,25:30], Inf)

norm((fd_datamat + datamat)[6:11,13:24], Inf)
norm((fd_datamat + datamat)[12:17,13:24], Inf)
norm((fd_datamat + datamat)[12:17,13:15], Inf)
norm((fd_datamat + datamat)[12:17,16:18], Inf)
norm((fd_datamat + datamat)[12:17,19:21], Inf)
norm((fd_datamat + datamat)[12:17,22:24], Inf)

norm((fd_datamat + datamat)[6:11,25:30], Inf)
norm((fd_datamat + datamat)[6:8,25:30], Inf)

norm((fd_datamat + datamat)[12:17,25:31], Inf)


(fd_datamat + datamat)[6:11,7:9]
(fd_datamat + datamat)[15:17,19:21]
fd_datamat[15:17,19:21]
-datamat[15:17,19:21]

ineqcs = collect(mech.contacts)
ineqcs[1].parentid
mech.bodies
ineqcs[2].parentid
ineqcs[3].parentid
ineqcs[4].parentid


collect(mech.joints)
collect(mech.bodies)
collect(mech.contacts)
5 + 2 * 6 + 4 * 8
12 * 2 + 6 + 1

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(1e10 * fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 8e-3
plot(Gray.(1e10 .* sensi))
plot(Gray.(fd_sensi))
