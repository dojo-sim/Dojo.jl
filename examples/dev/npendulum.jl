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

mech = getmechanism(:npendulum, Δt = 0.01, g = -9.81, Nlink = 2)
initialize!(mech, :npendulum, ϕ1 = 0.5)
initializeSimulation!(mech, true)
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]
body1.state.x1
body1.state.x2[1]
body1.state.q1
body1.state.q2[1]
body1.state.v15
body1.state.vsol[1]
body1.state.vsol[2]
body1.state.ϕ15
body1.state.ϕsol[1]
body1.state.ϕsol[2]

body2.state.x1
body2.state.x2[1]
body2.state.q1
body2.state.q2[1]
body2.state.v15
body2.state.vsol[1]
body2.state.vsol[2]
body2.state.ϕ15
body2.state.ϕsol[1]
body2.state.ϕsol[2]

storage = simulate!(mech, 1.11, record = true, solver = :mehrotra!)

visualize(mech, storage, vis = vis)

body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]
body1.state.x1
body1.state.x2[1]
body1.state.q1
body1.state.q2[1]
body1.state.v15
body1.state.vsol[1]
body1.state.vsol[2]
body1.state.ϕ15
body1.state.ϕsol[1]
body1.state.ϕsol[2]

body2.state.x1
body2.state.x2[1]
body2.state.q1
body2.state.q2[1]
body2.state.v15
body2.state.vsol[1]
body2.state.vsol[2]
body2.state.ϕ15
body2.state.ϕsol[1]
body2.state.ϕsol[2]

ex = [1.; 0; 0]
vert11 = 0*[0; 0; 1/2]
vert12 = -vert11
links = [Box(1, 1, 1, 1., color = RGBA(1., 0., 0.)) for i = 1:2]
ori = mech.origin
traxx = Prototype(:Revolute, ori, links[1], ex; p1 = vert12, p2 = vert11, spring = 0.0, damper = 0.0)[1][1]
jointb1 = EqualityConstraint(Prototype(:Revolute, ori, links[1], ex; p1 = vert12, p2 = vert11, spring = 0.0, damper = 0.0))
g(jointb1.constraints[1], SVector{3,Float64}(0,0,-0.5), one(UnitQuaternion))
g(traxx, SVector{3,Float64}(0,0,-0.5), one(UnitQuaternion))
g(tra1, SVector{3,Float64}(0,0,-0.5), one(UnitQuaternion))
g(tra2, SVector{3,Float64}(0,0,-0.5), one(UnitQuaternion))

tra1 = collect(mech.eqconstraints)[1].constraints[1]
rot1 = collect(mech.eqconstraints)[1].constraints[2]
tra2 = collect(mech.eqconstraints)[2].constraints[1]
rot2 = collect(mech.eqconstraints)[2].constraints[2]

tra1.vertices
tra2.vertices

################################################################################
# Differentiation
################################################################################

include(joinpath(module_dir(), "examples", "diff_tools.jl"))
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
@test norm(fd_datamat + datamat, Inf) < 1e-8
plot(Gray.(abs.(datamat)))
plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 1e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))
