# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())

using MeshCat

# Open visualizer
vis = Visualizer()
open(vis)

module_dir()
# Include new files
include(joinpath(module_dir(), "examples", "loader.jl"))


function getpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
        spring = 0.0, damper = 0.0, spring_offset = szeros(1), joint_limits = [-sones(1), sones(1)]) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1,
        joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset, rot_joint_limits = joint_limits))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end



mech = getmechanism(:pendulum, Δt = 0.01, g = 9.81)
initialize!(mech, :pendulum, ϕ1 = 0.1)
storage = simulate!(mech, 3.1, record = true, verbose = true)
visualize(mech, storage, vis=vis)




################################################################################
# Differentiation
################################################################################

# Set data
Nb = length(mech.bodies)
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

# # IFT
# datamat = full_data_matrix(mech)
setentries!(mech)
solmat = full_matrix(mech.system)
rank(solmat)
rank(solmat[1:9,1:9])
rank(solmat[1:9,1:9])
rank(solmat[10:15,10:15])
solmat[1:3,1:3]
solmat[1:9,1:9]

tra = eqc1.constraints[1]
rot = eqc1.constraints[2]
resetVars!(eqc1)
ηtra = eqc1.λsol[2][1:3]
ηrot = eqc1.λsol[2][3 .+ (1:6)]
∂g∂ʳself(mech, eqc1)
∂g∂ʳself(tra, ηtra)
∂g∂ʳself(rot, ηrot)

# sensi = - (solmat \ datamat)
# sensi2 = sensitivities(mech, sol, data)
#
# @test norm(sensi - sensi2, Inf) < 1.0e-8
# v0 = rand(13)
# @test norm(jvp(mech, sol, data, v0) - sensi * v0, Inf) < 1.0e-8
#
# # finite diff
# fd_datamat = finitediff_data_matrix(mech, data, sol, δ = 1e-5) * attjac
#
# @test norm(fd_datamat + datamat, Inf) < 1e-8
# plot(Gray.(abs.(datamat)))
# plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(1e11 .* solmat)))
plot(Gray.(abs.(1e9 .* fd_solmat)))

eqc1 = collect(mech.eqconstraints)[1]
λindex(eqc1, 1)
λindex(eqc1, 2)

∂g∂ʳpos()


fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr=1.0e-12, ϵb=1.0e-12) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 1e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

norm(fd_sensi - sensi, Inf) / norm(fd_sensi, Inf)
norm(fd_sensi - sensi, Inf)
