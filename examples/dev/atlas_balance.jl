# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))

# Open visualizer
vis = Visualizer()
open(vis)

# Build mechanism
mech = getmechanism(:atlas, Δt = 0.05, g = -0.4, cf = 0.8, contact = true)
initialize!(mech, :atlas, tran = [0,0,2.93], rot = [0.,0,0])
# storage = simulate!(mech, 2.0, record = true, solver = :mehrotra!)

# PD control law
nu = sum(getcontroldim.(mech.eqconstraints))
angles = [minimalCoordinates(mech, joint)[1] for joint in collect(mech.eqconstraints)[2:end]]
δangles = zeros(nu)
ind = 23
δangles[ind] += π/2
angles += δangles

function controller!(mechanism, k)
    for (i,joint) in enumerate(collect(mechanism.eqconstraints)[2:end])
        if getcontroldim(joint) == 1
            θ = minimalCoordinates(mechanism, joint)[1]
            dθ = minimalVelocities(mechanism, joint)[1]
            u = 1e-1 * (angles[i] - θ) + 5e-2 * (0 - dθ)
            u = clamp(u, -1.0, 1.0) * mechanism.Δt
            setForce!(mechanism, joint, SA[u])
        end
    end
    return
end

mech.bodies[33].state.vsol

forcedstorage = simulate!(mech, 0.50, controller!, record = true, solver = :mehrotra!)
visualize(mech, forcedstorage, vis = vis)

gains = zeros(30, 2)
gains[23,:] = [1e-1, 5e-2]

nams = [eqc.name for eqc in mech.eqconstraints]



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
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e10 * solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵ = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 8e-3
plot(Gray.(1e10 .* sensi))
plot(Gray.(fd_sensi))
