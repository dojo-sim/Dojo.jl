# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..", "..")
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


function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints)[1:end])
        nu = controldim(eqc)
        u = 0 * mechanism.Δt * ones(nu)
        setForce!(mechanism, eqc, u)
    end
    return
end


mech = getmechanism(:sphere, Δt = 0.05, g = -0*9.81);
initialize!(mech, :sphere, x=zeros(3), v=zeros(3), ω=zeros(3))
storage = simulate!(mech, 0.20, controller!, record=true, verbose=true)
visualize(mech, storage, vis=vis)



################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))1
# Set data
Nb = length(mech.bodies)
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

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

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr=1.0e-12, ϵb=1.0e-12) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 1e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

norm(fd_sensi - sensi, Inf) / norm(fd_sensi, Inf)
norm(fd_sensi - sensi, Inf)
