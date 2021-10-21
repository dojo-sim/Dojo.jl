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


# Build mechanism
include("mechanism_zoo.jl")


###############################################################################
# METHODS NEW EQUALITY CONSTRAINT
################################################################################

# t2r3
function TorqueRevolute(body1::AbstractBody{T}, body2, axis; p1 = szeros(T, 3), p2 = szeros(T, 3), qoffset = one(UnitQuaternion{T}), spring = zero(T), damper = zero(T)) where T
    return Translational3{T}(body1, body2; p1, p2, spring, damper),
    Rotational2{T}(body1, body2; axis, qoffset, spring, damper),
    Torque1{T}(body1, body2; axis=axis, qoffset=qoffset, spring=spring, damper=damper)
end

################################################################################
# DEVELOPMENT NEW EQUALITY CONSTRAINT
################################################################################
T = Float64

# Parameters
ex = [1.; 0; 0]
h = 1.
r = .05
vert11 = [0; 0; h/2]
vert12 = -vert11
Nlink = 3

# Links
origin = Origin{T}()
links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

# Constraints
spring0 = 1.0 * 1e1
damper0 = 1.0 * 1e2
spring1 = 1.0 * 1e1
damper1 = 1.0 * 1e2
jointb1 = EqualityConstraint(Revolute(origin, links[1], ex; spring=spring0, damper=damper0, p2 = vert11))
if Nlink > 1
    eqcs = [
        jointb1;
        [EqualityConstraint(TorqueRevolute(links[i - 1], links[i], ex; spring=spring1, damper=damper1, p1=vert12, p2=vert11)) for i = 2:Nlink]
        ]
else
    eqcs = [jointb1]
end
# jointb1 = EqualityConstraint(TorqueRevolute(origin, links[1], ex; spring=spring0, damper=damper0, p2 = vert11))
# if Nlink > 1
#     eqcs = [
#         jointb1;
#         [EqualityConstraint(Revolute(links[i - 1], links[i], ex; spring=spring1, damper=damper1, p1=vert12, p2=vert11)) for i = 2:Nlink]
#         ]
# else
#     eqcs = [jointb1]
# end
mech = Mechanism(origin, links, eqcs, g = -9.81, Δt = 0.01)

initialize!(mech, :npendulum)
storage = simulate!(mech, 1.0, record = true, solver = :mehrotra!)

# visualize(mech, storage, vis = vis)
# tor = jointb1.constraints[3]
# xb = mech.bodies[2].state.xc
# qb = mech.bodies[2].state.qc
# ωb = mech.bodies[2].state.ωc
# Δt = mech.Δt
# fg = q -> ∂g∂ʳposb(tor, xb, UnitQuaternion(q...))[1:3, 4:6]# #* LVᵀmat(UnitQuaternion(q...))
# fg(qb)
# ForwardDiff.jacobian(fg, [qb.w; qb.x; qb.y; qb.z])

################################################################################
# Differentiation
################################################################################
include(joinpath(module_dir(), "examples", "dev", "diff_tools.jl"))
# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)

# eqc1 = collect(mech.eqconstraints)[1]
# eqc2 = collect(mech.eqconstraints)[1]
# tra1 = eqc1.constraints[1]
# rot1 = eqc1.constraints[2]
# torque1 = eqc1.constraints[3]
# tra2 = eqc2.constraints[1]
# rot2 = eqc2.constraints[2]
# torque2 = eqc2.constraints[3]
# constraintmat(tra1)
# constraintmat(rot1)
# constraintmat(torque1)
# constraintmat(tra2)
# constraintmat(rot2)
# constraintmat(torque2)

# IFT
setentries!(mech)
datamat = full_data_matrix(mech)
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)
@show norm(full_vector(mech.system), Inf)

# finite diff
fd_datamat = finitediff_data_matrix(mech, data, sol) * attjac
-datamat
@test norm(fd_datamat + datamat, Inf) < 1e-7
# plot(Gray.(abs.(datamat)))
# plot(Gray.(abs.(fd_datamat)))
norm((fd_datamat + datamat)[1:6, 1:13], Inf)
norm((fd_datamat + datamat)[7:12, 1:13], Inf)
norm((fd_datamat + datamat)[7:12, 1:3], Inf)
norm((fd_datamat + datamat)[7:12, 4:6], Inf)
norm((fd_datamat + datamat)[7:12, 7:10], Inf)
norm((fd_datamat + datamat)[7:9, 7:10], Inf)
norm((fd_datamat + datamat)[10:12, 7:10], Inf)
norm((fd_datamat + datamat)[7:12, 11:13], Inf)

norm((fd_datamat + datamat)[1:12, 1:26], Inf)
norm((fd_datamat + datamat)[13:18, 1:13], Inf)
norm((fd_datamat + datamat)[13:18, 14:26], Inf)
datamat[13:18, 24:26]
-fd_datamat[13:18, 24:26]

norm((fd_datamat + datamat)[19:24, 1:26], Inf)
norm((fd_datamat + datamat)[19:24, 1:13], Inf)
norm((fd_datamat + datamat)[19:24, 14:26], Inf)
datamat[19:24, 24:26]
-fd_datamat[19:24, 24:26]


(fd_datamat + datamat)[10:12, 7:10]
fd_datamat[10:12, 7:10]
datamat[10:12, 7:10]

fd_solmat = finitediff_sol_matrix(mech, data, sol)
@test norm(fd_solmat + solmat, Inf) < 1e-7
plot(Gray.(abs.(solmat)))
plot(Gray.(abs.(fd_solmat)))

fd_sensi = finitediff_sensitivity(mech, data) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 5e-3
# plot(Gray.(sensi))
# plot(Gray.(fd_sensi))
# norm(fd_sensi - sensi, Inf)
# norm(fd_sensi, Inf)
# norm(fd_sensi - sensi) / norm(fd_sensi)
