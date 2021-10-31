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

mech = getmechanism(:atlas, Δt = 0.01, g = -2.0, cf = 0.8, contact = true, spring = 1000.0, damper = 50.0)
initialize!(mech, :atlas, tran = [0,0,1.0], rot = [0.0,0,0])

function controller!(mechanism, k)
    for (i,eqc) in enumerate(collect(mechanism.eqconstraints)[2:end])
        pbody = getbody(mech, eqc.parentid)
        minJ = minimum(diag(pbody.J))
        for (i,joint) in enumerate(eqc.constraints)
            cbody = getbody(mech, eqc.childids[i])
            minJ = min(minJ, minimum(diag(cbody.J)))
        end
        nu = getcontroldim(eqc)
        u = 1 * minJ * (rand(nu) .- 0.2) * Δt_ * 0.0
        setForce!(mechanism, eqc, SVector{nu}(u))
    end
    return
end

@profiler storage = simulate!(mech, 1.50, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

[eqc.isspring for eqc in collect(mech.eqconstraints)]

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
