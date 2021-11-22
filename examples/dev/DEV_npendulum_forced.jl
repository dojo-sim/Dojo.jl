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

# Build mechanism
mech = getmechanism(:npendulum, Δt = 0.05, g = -9.81, Nlink = 5)
initialize!(mech, :npendulum, ϕ1 = 0.1*pi)

function cont!(mechanism, k; u = 1.0)
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        nu = controldim(eqc, ignore_floating_base = false)
        su = mechanism.Δt * u * sones(nu)
        setForce!(mechanism, eqc, su)
    end
    return
end

storage = simulate!(mech, 10.0, cont!, record = true, solver = :mehrotra!)
visualize(mech, storage, vis = vis)
plot([q.x for q in storage.q[1]])
control_datamat(mech)

################################################################################
# Differentiation
################################################################################

# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)

# IFT
setentries!(mech)
datamat = full_data_matrix(deepcopy(mech))
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(deepcopy(mech), data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-7
plot(Gray.(abs.(1e8 * datamat)))
plot(Gray.(abs.(1e8*fd_datamat)))

norm((datamat + fd_datamat)[1:3, 1:15], Inf)
norm((datamat + fd_datamat)[4:9, 1:12], Inf)
norm((datamat + fd_datamat)[4:9, 13:15], Inf)

# norm((datamat + fd_datamat)[6:10, 13], Inf)

datamat[4:9, 13:15]
-fd_datamat[4:9, 13:15]

norm((datamat + fd_datamat)[11:16, 1:26], Inf)
norm((datamat + fd_datamat)[17:22, 1:26], Inf)

norm((datamat + fd_datamat)[6:11, 1:13], Inf)

norm((datamat + fd_datamat)[6:11, 1:3], Inf)
norm((datamat + fd_datamat)[6:11, 4:6], Inf)
norm((datamat + fd_datamat)[6:11, 7:9], Inf)
norm((datamat + fd_datamat)[6:11, 10:12], Inf)
norm((datamat + fd_datamat)[6:11, 13:13], Inf)

datamat[6:11, 1:3]
-fd_datamat[6:11, 1:3]

joint = mech.eqconstraints[1].constraints[1]
xb = mech.bodies[2].state.x2[1]
qb = mech.bodies[2].state.q2[1]

length(mech.eqconstraints[1])
length(joint)
fw = w -> transpose(∂g∂ʳposb(joint, w, qb)[1:3, 4:6]) * constraintmat(joint)' * mech.eqconstraints[1].λsol[2][1:3]
ForwardDiff.jacobian(fw, xb)

using FiniteDiff
fz = z -> -transpose(∂g∂ʳposb(joint, xb, UnitQuaternion(z..., false))[1:3, 1:3]) * constraintmat(joint)' * mech.eqconstraints[1].λsol[2][1:3]
FiniteDiff.finite_difference_jacobian(fz, [qb.w, qb.x, qb.y, qb.z]) * LVᵀmat(qb)
ForwardDiff.jacobian(fz, [qb.w, qb.x, qb.y, qb.z]) #* LVᵀmat(qb)
fdjac(fz, [qb.w, qb.x, qb.y, qb.z]) #* LVᵀmat(qb)

norm((datamat + fd_datamat)[1:10, 1:26], Inf)
norm((datamat + fd_datamat)[11:16, 1:6], Inf)
norm((datamat + fd_datamat)[11:16, 4:6], Inf)
norm((datamat + fd_datamat)[11:16, 7:9], Inf)
norm((datamat + fd_datamat)[11:16, 10:12], Inf)

norm((datamat + fd_datamat)[11:16, 13:15], Inf)
norm((datamat + fd_datamat)[11:16, 16:18], Inf)
norm((datamat + fd_datamat)[11:16, 19:21], Inf)
norm((datamat + fd_datamat)[11:16, 22:24], Inf)
norm((datamat + fd_datamat)[11:16, 25:26], Inf)

norm((datamat + fd_datamat)[17:22, 1:6], Inf)
norm((datamat + fd_datamat)[17:22, 4:6], Inf)
norm((datamat + fd_datamat)[17:22, 7:9], Inf)
norm((datamat + fd_datamat)[17:22, 10:12], Inf)

norm((datamat + fd_datamat)[17:22, 13:15], Inf)
norm((datamat + fd_datamat)[17:22, 16:18], Inf)
norm((datamat + fd_datamat)[17:22, 19:21], Inf)
norm((datamat + fd_datamat)[17:22, 22:24], Inf)
norm((datamat + fd_datamat)[17:22, 25:26], Inf)

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
plot(Gray.(10e8 * abs.(solmat)))
plot(Gray.(10e8 * abs.(fd_solmat)))
norm(fd_solmat + solmat, Inf)
norm((fd_solmat + solmat)[1:10, 1:10], Inf)
norm((fd_solmat + solmat)[1:10, 11:22], Inf)
norm((fd_solmat + solmat)[11:22, 1:10], Inf)
norm((fd_solmat + solmat)[11:22, 11:22], Inf)
fd_solmat[11:16, 17:22]
-solmat[11:16, 17:22]

fd_solmat[17:22, 11:16]
-solmat[17:22, 11:16]



norm((fd_solmat + solmat)[11:16, 11:16], Inf)
norm((fd_solmat + solmat)[11:16, 17:22], Inf)
norm((fd_solmat + solmat)[17:22, 11:16], Inf)
norm((fd_solmat + solmat)[17:22, 17:22], Inf)

(fd_solmat + solmat)[11:16, 11:16]
(fd_solmat + solmat)[11:16, 11:16][4:6,4:6]
(fd_solmat + solmat)[17:22, 11:16]
(fd_solmat + solmat)[17:22, 17:22]

fd_solmat[11:16, 11:16]
fd_solmat[11:16, 11:16][4:6, 4:6]
fd_solmat[11:16, 17:22]
fd_solmat[17:22, 11:16]
fd_solmat[17:22, 17:22]

solmat[11:16, 11:16]
solmat[11:16, 11:16][4:6, 4:6]
solmat[11:16, 17:22]
solmat[17:22, 11:16]
solmat[17:22, 17:22]

norm(solmat, Inf)







################################################################################
# Finite Diff
################################################################################
include(joinpath(@__DIR__, "finite_diff.jl"))
mech.Δt
Δt = 0.01
rot1 = mech.eqconstraints[1].constraints[2]
rot2 = mech.eqconstraints[2].constraints[2]
origin = mech.origin
body1 = collect(mech.bodies)[1]
body2 = collect(mech.bodies)[2]

jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, springforcea, ∂springforcea∂vela, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, damperforcea, ∂damperforcea∂vela, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, springforcea, ∂springforcea∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, damperforcea, ∂damperforcea∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, springforceb, ∂springforceb∂vela, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot2, body1, body2, Δt, damperforceb, ∂damperforceb∂vela, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_vel(rot1, origin, body1, Δt, springforceb, ∂springforceb∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0
jac1
jac0, jac1 = finitediff_vel(rot1, origin, body1, Δt, damperforceb, ∂damperforceb∂velb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8


jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, springforcea, ∂springforcea∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, damperforcea, ∂damperforcea∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, springforcea, ∂springforcea∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, damperforcea, ∂damperforcea∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, springforceb, ∂springforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, damperforceb, ∂damperforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, springforceb, ∂springforceb∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot2, body1, body2, Δt, damperforceb, ∂damperforceb∂posa, diff_body = :parent)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot1, origin, body1, Δt, springforceb, ∂springforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8
jac0, jac1 = finitediff_pos(rot1, origin, body1, Δt, damperforceb, ∂damperforceb∂posb, diff_body = :child)
@test norm(jac0 - jac1, Inf) < 1e-8







# solmat[1:5, 1:5]
# solmat[1:5, 6:11]
# solmat[6:11, 1:5]
# solmat[6:11, 6:11]
# solmat[9:11, 9:11]
#
#
#
# fd_solmat[1:5, 1:5]
# fd_solmat[1:5, 6:11]
# fd_solmat[6:11, 1:5]
# fd_solmat[6:11, 6:11]
# fd_solmat[9:11, 9:11]
#
#
#
# (solmat + fd_solmat)[1:5, 1:5]
# (solmat + fd_solmat)[1:5, 6:11]
# (solmat + fd_solmat)[6:11, 1:5]
# (solmat + fd_solmat)[6:11, 6:11]
# (solmat + fd_solmat)[9:11, 9:11]

fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 3e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

# diagonal∂damper∂ʳvel(mech.eqconstraints[1],
# offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# diagonal∂damper∂ʳvel(mech, mech.eqconstraints[1], mech.bodies[2])
# offdiagonal∂damper∂ʳvel(mech.eqconstraints[1].constraints[1], mech.origin, mech.bodies[2], mech.bodies[2].id, mech.Δt)
# offdiagonal∂damper∂ʳvel(mech.eqconstraints[1].constraints[2], mech.origin, mech.bodies[2], mech.bodies[2].id, mech.Δt)

# mech.bodies

# ################################################################################
# # Damper Jacobian
# ################################################################################

# include("fd_tools.jl")

# j0 = mech.eqconstraints[1]
# jt0 = j0.constraints[1]
# jr0 = j0.constraints[2]
# origin0 = mech.origin
# bodya0 = mech.bodies[3]
# bodyb0 = mech.bodies[4]
# childida0 = 3
# childidb0 = 4
# Δt0 = mech.Δt
# damperforcea(jt0, bodya0, bodyb0, childidb0, Δt0)
# damperforceb(jt0, bodya0, bodyb0, childidb0, Δt0)
# damperforceb(jr0, origin0, bodya0, childida0, Δt0)

# x2a0, q2a0 = posargs3(bodya0.state, Δt0)
# x2b0, q2b0 = posargs3(bodyb0.state, Δt0)
# x1a0, v1a0, q1a0, ω1a0 = fullargssol(bodya0.state)
# x1b0, v1b0, q1b0, ω1b0 = fullargssol(bodyb0.state)

# Random.seed!(100)
# x2a0 = rand(3)
# q2a0 = UnitQuaternion(rand(4)...)
# x2b0 = rand(3)
# q2b0 = UnitQuaternion(rand(4)...)
# x1a0 = rand(3)
# v1a0 = rand(3)
# q1a0 = UnitQuaternion(rand(4)...)
# ω1a0 = rand(3)
# x1b0 = rand(3)
# v1b0 = rand(3)
# q1b0 = UnitQuaternion(rand(4)...)
# ω1b0 = rand(3)



# ################################################################################
# # Damper translation
# ################################################################################
# jt0 = FixedOrientation(bodya0, bodyb0; spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)

# # jt0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)
# #
# # jt0 = Prismatic(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)
# #
# # jt0 = Fixed(bodya0, bodyb0)[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)

# damperforcea(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# damperforceb(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# damperforceb(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# Dtra1 = diagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Dtra2 = offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Dtra3 = offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# fd_Dtra1 = fd_diagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Dtra2 = fd_offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Dtra3 = fd_offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# norm(Dtra1 - fd_Dtra1)
# norm(Dtra2 - fd_Dtra2)
# norm(Dtra3 - fd_Dtra3)



# ################################################################################
# # Damper rotation
# ################################################################################
# jr0 = Spherical(bodya0, bodyb0, spring = zeros(3), damper = zeros(3))[2][1]
# jr0.spring = 1e1 .* rand(3)
# jr0.damper = 1e1 .* rand(3)

# # # jr0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # # jr0.spring = 1e1 .* rand(3)
# # # jr0.damper = 1e1 .* rand(3)
# #
# # jr0 = Revolute(bodya0, bodyb0, rand(3))[2][1]
# # jr0.spring = 1e1 .* rand(3)
# # jr0.damper = 1e1 .* rand(3)
# #
# # jr0 = Fixed(bodya0, bodyb0)[2][1]
# # jr0.spring = 1e1 .* rand(3)
# # jr0.damper = 1e1 .* rand(3)

# damperforcea(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# damperforceb(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# damperforceb(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# Drot1 = diagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Drot2 = offdiagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Drot3 = offdiagonal∂damper∂ʳvel(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# fd_Drot1 = fd_diagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Drot2 = fd_offdiagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Drot3 = fd_offdiagonal∂damper∂ʳvel(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# norm(Drot1 - fd_Drot1)
# norm(Drot2 - fd_Drot2)
# norm(Drot3 - fd_Drot3)

# function der1(ω1a, q2a, ω1b, q2b)
#     invqbqa = q2b\q2a
#     A = nullspacemat(jr0)
#     AᵀA = zerodimstaticadjoint(A) * A
#     return 2*VLmat(invqbqa)*RVᵀmat(invqbqa)* AᵀA * Diagonal(jr0.damper) * AᵀA
# end

# function der3(ω1a, q2a, ω1b, q2b)
#     A = I(3)
#     Aᵀ = A'
#     C = -2 * Aᵀ * A * Diagonal(jr0.damper) * Aᵀ * A
#     Δq = q2a \ q2b
#     Δqbar = q2b \ q2a
#     dF1 = C * VRᵀmat(Δq) * LVᵀmat(Δq)
#     dF2 = VRᵀmat(Δqbar) * LVᵀmat(Δqbar) * dF1
#     return dF2
# end

# function der4(ω1a, q2a, ω1b, q2b)
#     A = I(3)
#     Aᵀ = A'
#     function f(ω1b)
#         q2a_ = UnitQuaternion(q2a.w, q2a.x, q2a.y, q2a.z, false)
#         q2b_ = UnitQuaternion(q2b.w, q2b.x, q2b.y, q2b.z, false)
#         velocity = A * (vrotate(ω1b,q2a_\q2b_) - ω1a) # in body1's frame
#         force = -2 * Aᵀ * A * Diagonal(jr0.damper) * Aᵀ * velocity
#         force = vrotate(force, q2b_ \ q2a_) # in body2's frame
#         return force
#     end
#     ForwardDiff.jacobian(f, ω1b)
# end

# ω1a = rand(3)
# q2a = UnitQuaternion(rand(4)...)
# ω1b = rand(3)
# q2b = UnitQuaternion(rand(4)...)
# d1 = der1(ω1a, q2a, ω1b, q2b)
# d3 = der3(ω1a, q2a, ω1b, q2b)
# d4 = der4(ω1a, q2a, ω1b, q2b)
# norm(d4 - d3)
# norm(d4 - d1)

# ################################################################################
# # Spring translation
# ################################################################################
# jt0 = FixedOrientation(bodya0, bodyb0; spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)

# # jt0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)
# #
# # jt0 = Prismatic(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)
# #
# # jt0 = Fixed(bodya0, bodyb0)[1][1]
# # jt0.spring = 1e1 .* rand(3)
# # jt0.damper = 1e1 .* rand(3)

# springforcea(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# springforceb(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# springforceb(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# Dspr1 = diagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Dspr2 = offdiagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# Dspr3 = offdiagonal∂spring∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# fd_Dspr1 = fd_diagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Dspr2 = fd_offdiagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
# fd_Dspr3 = fd_offdiagonal∂spring∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

# norm(Dspr1 - fd_Dspr1)
# norm(Dspr2 - fd_Dspr2)
# (Dspr2 - fd_Dspr2)[1:3,1:3]
# (Dspr2 - fd_Dspr2)[1:3,1:3]
# norm(Dspr3 - fd_Dspr3)

# Z = szeros(3, 3)
# Z1 = sones(3, 3)
# [[Z1; Z] [Z1; Z]]
