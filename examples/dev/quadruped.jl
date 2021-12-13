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


mech = getmechanism(:quadruped, Δt = 0.01, g = -9.81, cf = 0.8, contact = true, damper = 1.0, spring = 30.0)
initialize!(mech, :quadruped, tran = [0,0,0.56], rot = [0.10,0.05,0.03], θ = 0.95)
@elapsed storage = simulate!(mech, 0.5, record = true, verbose = false)
visualize(mech, storage, vis = vis)

Δt_ = 0.001
begin
    N = 100
    steps = Base.OneTo(N)
    storage = Storage(steps,length(bodies))

    traj21 = [initangle*(cos(i*0.001*2*pi)*0.1+0.9) for i=1:N]
    traj21 = [[initangle for i=1:500];traj21]
    traj21 = SA[traj21...]
    traj31 = [-2*initangle-sin(i*0.001*2*pi)*0.1 for i=1:N]
    traj31 = [[-2*initangle for i=1:500];traj31]
    traj31 = SA[traj31...]

    traj22 = [initangle*(cos(i*0.001*2*pi+pi)*0.1+0.9) for i=1:N]
    traj22 = [[initangle*0.9 for i=1:500];traj22]
    traj22 = SA[traj22...]
    traj32 = [-2*initangle-sin(i*0.001*2*pi+pi)*0.1 for i=1:N]
    traj32 = [[-2*initangle for i=1:500];traj32]
    traj32 = SA[traj32...]



    function singleleg(mechanism, leg, angles)
        j1 = geteqconstraint(mechanism, leg*"_hip_joint")
        j2 = geteqconstraint(mechanism, leg*"_thigh_joint")
        j3 = geteqconstraint(mechanism, leg*"_calf_joint")

        θ1 = minimalCoordinates(mechanism, j1)[1]
        θ2 = minimalCoordinates(mechanism, j2)[1]
        θ3 = minimalCoordinates(mechanism, j3)[1]
        dθ1 = minimalVelocities(mechanism, j1)[1]
        dθ2 = minimalVelocities(mechanism, j2)[1]
        dθ3 = minimalVelocities(mechanism, j3)[1]

        u1 = (100.0*(angles[1]-θ1) + 5.0*(0-dθ1)) * Δt_ #* 0.17
        u2 = (80.0*(angles[2]-θ2) + 4.0*(0-dθ2)) * Δt_ #* 0.17
        u3 = (60.0*(angles[3]-θ3) + 3.0*(0-dθ3)) * Δt_ #* 0.17

        setForce!(mechanism, j1, SA[u1])
        setForce!(mechanism, j2, SA[u2])
        setForce!(mechanism, j3, SA[u3])
    end

    function controller!(mechanism, k)
        singleleg(mechanism, "FR", SA[0.0;traj21[k];traj31[k]])
        singleleg(mechanism, "FL", SA[0.0;traj22[k];traj32[k]])
        singleleg(mechanism, "RR", SA[0.0;traj22[k];traj32[k]])
        singleleg(mechanism, "RL", SA[0.0;traj21[k];traj31[k]])
    end
end



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
