################################################################################
# Development
################################################################################
using ConstrainedDynamics
using ConstrainedDynamicsVis

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
using StaticArrays
using LinearAlgebra
using Test


# Include dev files
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))


Δt = 0.01
g0 = -9.81

# Parameters
joint_axis = [1.0;0.0;0.0]
length1 = 0.5
width, depth = 0.5, 0.5

origin = Origin{Float64}()
link1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
joint0to1 = EqualityConstraint(Floating(origin, link1))
links = [link1]
eqcs = [joint0to1]

# Corner vectors
corners = [
    # [[0 ; 0; 0.]]
    [[length1 / 2;length1 / 2;-length1 / 2]]
    [[length1 / 2;-length1 / 2;-length1 / 2]]
    # [[-length1 / 2;length1 / 2;-length1 / 2]]
    # [[-length1 / 2;-length1 / 2;-length1 / 2]]
    # [[length1 / 2;length1 / 2;length1 / 2]]
    # [[length1 / 2;-length1 / 2;length1 / 2]]
    # [[-length1 / 2;length1 / 2;length1 / 2]]
    # [[-length1 / 2;-length1 / 2;length1 / 2]]
]
n = length(corners)
normal = [[0;0;1.0] for i = 1:n]
cf = 0.2 * ones(n)

contbounds = [ContactBound(links[1], normal[i], cf[i], p = corners[i]) for i=1:length(corners)]
contineqcs = [InequalityConstraint((contbounds[i], links[1].id, nothing)) for i=1:length(corners)]

# fricsandineqs = [Friction(link1, normal[i], cf[i]; p = corners[i]) for i=1:n]
# frics = getindex.(fricsandineqs,1)
# ineqcs = vcat(getindex.(fricsandineqs,2)...)
# imps = ineqcs[1:3:end]
# mech = Mechanism(origin, links, eqcs, imps, g = g0, Δt = Δt)
mech = Mechanism(origin, links, eqcs, contineqcs, g = g0, Δt = Δt)



Random.seed!(100)
ω = 1.0 * (rand(3) .- 0.5) * 1
x = [0, 0, 1.0]
v = [1, 0.3, 0.2]
initialize!(mech, :dice, x = x, v = v, ω = ω)
storage = simulate!(mech, 0.5, record = true, solver = :mehrotra!)

plot(hcat(Vector.(storage.x[1])...)')
# visualize(mech, storage)



include(joinpath(module_dir(), "examples", "dev", "diff_tools_control_contact.jl"))
# Set data
Nb = length(mech.bodies)
data = getdata(mech)
setdata!(mech, data)
mehrotra!(mech, opts = InteriorPointOptions(rtol = 1e-12, btol = 1e-12))
sol = getsolution(mech)
attjac = attitudejacobian(data, Nb)

# IFT
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)


# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

# fd_datamat[1:6,1:6]
# fd_datamat[1:6,7:12]
# fd_datamat[7:14,1:6]
# fd_datamat[7:14,7:12]
# fd_datamat[11:14,1:6]
# fd_datamat[11:14,7:12]
#
# datamat[1:6,1:6]
# datamat[1:6,7:12]
# datamat[7:14,1:6]
# datamat[7:14,7:12]
# datamat[11:14,1:6]
# datamat[11:14,7:12]
#
# (fd_datamat + datamat)[1:6,1:6]
# (fd_datamat + datamat)[1:6,7:12]
# (fd_datamat + datamat)[4:6,7:12]
# (fd_datamat + datamat)[4:6,7:9]
# (fd_datamat + datamat)[4:6,10:12]
# (fd_datamat + datamat)[7:14,1:6]
# (fd_datamat + datamat)[7:14,7:12]
# (fd_datamat + datamat)[11:14,1:3]
# (fd_datamat + datamat)[11:14,7:9]
# norm(fd_datamat + datamat, Inf)

fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat) < 1e-7
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
fd2_sensi = fd_solmat \ fd_datamat
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 6e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

norm(fd_sensi - sensi)
norm(fd_sensi - sensi) / norm(fd_sensi)
