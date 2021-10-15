# Utils
function module_dir()
    return joinpath(@__DIR__, "..", "..")
end

# Activate package
using Pkg
Pkg.activate(module_dir())
include(joinpath(module_dir(), "examples", "dev", "loader.jl"))

# Open visualizer
# vis = Visualizer()
# open(vis)

# Build mechanism
mech = getmechanism(:npendulum, Δt = 0.01, g = -9.81, Nlink = 2)
initialize!(mech, :npendulum, ϕ1 = 1.3)

for (i,joint) in enumerate(mech.eqconstraints)
    if i ∈ (1,2)
        jt = joint.constraints[1]
        jr = joint.constraints[2]
        joint.isdamper = true #false
        joint.isspring = true #false

        jt.spring = 1/i * 0.0 * 1e-0# 1e4
        jt.damper = 1/i * 3.3 * 1e+3# 1e4
        jr.spring = 1/i * 0.0 * 1e-0# 1e4
        jr.damper = 1/i * 3.3 * 1e+3# 1e4

        mech.eqconstraints[1].isspring
        mech.eqconstraints[1].isdamper
        mech.eqconstraints[1].constraints[2].damper
    end
end

storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)
# visstorage = simulate!(mech, 4.0, record = true, solver = :mehrotra!)
# plot(hcat(Vector.(storage.x[1])...)')
# plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
# plot(hcat(Vector.(storage.v[1])...)')
# plot(hcat(Vector.(storage.ω[1])...)')

# visualize(mech, visstorage, vis = vis)


################################################################################
# Differentiation
################################################################################

# Set data
data = getdata(mech)
setdata!(mech, data)
sol = getsolution(mech)
Nb = length(collect(mech.bodies))
attjac = attitudejacobian(data, Nb)


setentries!(mech)
# IFT
datamat = full_data_matrix(deepcopy(mech))
solmat = full_matrix(mech.system)
sensi = - (solmat \ datamat)

# finite diff
fd_datamat = finitediff_data_matrix(deepcopy(mech), data, sol, δ = 1e-5) * attjac
@test norm(fd_datamat + datamat, Inf) < 1e-8
# plot(Gray.(abs.(datamat)))
# plot(Gray.(abs.(fd_datamat)))

fd_solmat = finitediff_sol_matrix(mech, data, sol, δ = 1e-5)
@test norm(fd_solmat + solmat, Inf) < 1e-8
# plot(Gray.(abs.(solmat)))
# plot(Gray.(abs.(fd_solmat)))
norm(fd_solmat + solmat, Inf)

norm((fd_solmat + solmat)[1:10, 1:10], Inf)
norm((fd_solmat + solmat)[1:10, 11:22], Inf)

norm((fd_solmat + solmat)[11:16, 11:16], Inf)
norm((fd_solmat + solmat)[11:16, 17:22], Inf)

norm((fd_solmat + solmat)[17:22, 11:16], Inf)
norm((fd_solmat + solmat)[17:22, 17:22], Inf)

norm((fd_solmat + solmat)[11:22, 1:10], Inf)
norm((fd_solmat + solmat)[11:22, 11:22], Inf)


mech.eqconstraints[1].constraints[1] 
mech.eqconstraints[1].constraints[2] 


solmat[11:16, 17:22]
fd_solmat[11:16, 17:22]

solmat[17:22, 11:16]
fd_solmat[17:22, 11:16]

∂gab∂ʳba(mech, mech.bodies[3], mech.bodies[4])[2]

offdiagonal∂damper∂ʳvel(mech.eqconstraints[1].constraints[2], mech.bodies[3].state.xsol[1], mech.bodies[3].state.qsol[1], mech.bodies[4].state.xsol[1], mech.bodies[4].state.qsol[1])
offdiagonal∂damper∂ʳvel(mech, mech.eqconstraints[1], mech.bodies[4], mech.bodies[3])

solmat[1:5, 6:11]
solmat[6:11, 1:5]
solmat[6:11, 6:11]
solmat[9:11, 9:11]

setentries!(mech) 

fd_solmat[1:5, 1:5]
fd_solmat[1:5, 6:11]
fd_solmat[6:11, 1:5]
fd_solmat[6:11, 6:11]
fd_solmat[9:11, 9:11]



(solmat + fd_solmat)[1:5, 1:5]
(solmat + fd_solmat)[1:5, 6:11]
(solmat + fd_solmat)[6:11, 1:5]
(solmat + fd_solmat)[6:11, 6:11]
(solmat + fd_solmat)[9:11, 9:11]



fd_sensi = finitediff_sensitivity(mech, data, δ = 1e-5, ϵr = 1e-14, ϵb = 1e-14) * attjac
@test norm(fd_sensi - sensi) / norm(fd_sensi) < 3e-3
plot(Gray.(sensi))
plot(Gray.(fd_sensi))

diagonal∂damper∂ʳvel(mech.eqconstraints[1],
offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)
diagonal∂damper∂ʳvel(mech, mech.eqconstraints[1], mech.bodies[2])
offdiagonal∂damper∂ʳvel(mech.eqconstraints[1].constraints[1], mech.origin, mech.bodies[2], mech.bodies[2].id, mech.Δt)
offdiagonal∂damper∂ʳvel(mech.eqconstraints[1].constraints[2], mech.origin, mech.bodies[2], mech.bodies[2].id, mech.Δt)

mech.bodies

################################################################################
# Damper Jacobian
################################################################################

include("fd_tools.jl")

j0 = mech.eqconstraints[1]
jt0 = j0.constraints[1]
jr0 = j0.constraints[2]
origin0 = mech.origin
bodya0 = mech.bodies[3]
bodyb0 = mech.bodies[4]
childida0 = 3
childidb0 = 4
Δt0 = mech.Δt
damperforcea(jt0, bodya0, bodyb0, childidb0, Δt0)
damperforceb(jt0, bodya0, bodyb0, childidb0, Δt0)
damperforceb(jr0, origin0, bodya0, childida0, Δt0)

x2a0, q2a0 = posargsnext(bodya0.state, Δt0)
x2b0, q2b0 = posargsnext(bodyb0.state, Δt0)
x1a0, v1a0, q1a0, ω1a0 = fullargssol(bodya0.state)
x1b0, v1b0, q1b0, ω1b0 = fullargssol(bodyb0.state)

Random.seed!(100)
x2a0 = rand(3)
q2a0 = UnitQuaternion(rand(4)...)
x2b0 = rand(3)
q2b0 = UnitQuaternion(rand(4)...)
x1a0 = rand(3)
v1a0 = rand(3)
q1a0 = UnitQuaternion(rand(4)...)
ω1a0 = rand(3)
x1b0 = rand(3)
v1b0 = rand(3)
q1b0 = UnitQuaternion(rand(4)...)
ω1b0 = rand(3)



################################################################################
# Damper translation
################################################################################
jt0 = FixedOrientation(bodya0, bodyb0; spring = zeros(3), damper = zeros(3))[1][1]
jt0.spring = 1e1 .* rand(3)
jt0.damper = 1e1 .* rand(3)

# jt0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)
#
# jt0 = Prismatic(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)
#
# jt0 = Fixed(bodya0, bodyb0)[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)

damperforcea(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

Dtra1 = diagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Dtra2 = offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Dtra3 = offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

fd_Dtra1 = fd_diagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dtra2 = fd_offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dtra3 = fd_offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

norm(Dtra1 - fd_Dtra1)
norm(Dtra2 - fd_Dtra2)
norm(Dtra3 - fd_Dtra3)



################################################################################
# Damper rotation
################################################################################
jr0 = Spherical(bodya0, bodyb0, spring = zeros(3), damper = zeros(3))[2][1]
jr0.spring = 1e1 .* rand(3)
jr0.damper = 1e1 .* rand(3)

# # jr0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# # jr0.spring = 1e1 .* rand(3)
# # jr0.damper = 1e1 .* rand(3)
#
# jr0 = Revolute(bodya0, bodyb0, rand(3))[2][1]
# jr0.spring = 1e1 .* rand(3)
# jr0.damper = 1e1 .* rand(3)
#
# jr0 = Fixed(bodya0, bodyb0)[2][1]
# jr0.spring = 1e1 .* rand(3)
# jr0.damper = 1e1 .* rand(3)

damperforcea(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
damperforceb(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

Drot1 = diagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Drot2 = offdiagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Drot3 = offdiagonal∂damper∂ʳvel(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

fd_Drot1 = fd_diagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Drot2 = fd_offdiagonal∂damper∂ʳvel(jr0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Drot3 = fd_offdiagonal∂damper∂ʳvel(jr0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

norm(Drot1 - fd_Drot1)
norm(Drot2 - fd_Drot2)
norm(Drot3 - fd_Drot3)

function der1(ω1a, q2a, ω1b, q2b)
    invqbqa = q2b\q2a
    A = nullspacemat(jr0)
    AᵀA = zerodimstaticadjoint(A) * A
    return 2*VLmat(invqbqa)*RVᵀmat(invqbqa)* AᵀA * Diagonal(jr0.damper) * AᵀA
end

function der3(ω1a, q2a, ω1b, q2b)
    A = I(3)
    Aᵀ = A'
    C = -2 * Aᵀ * A * Diagonal(jr0.damper) * Aᵀ * A
    Δq = q2a \ q2b
    Δqbar = q2b \ q2a
    dF1 = C * VRᵀmat(Δq) * LVᵀmat(Δq)
    dF2 = VRᵀmat(Δqbar) * LVᵀmat(Δqbar) * dF1
    return dF2
end

function der4(ω1a, q2a, ω1b, q2b)
    A = I(3)
    Aᵀ = A'
    function f(ω1b)
        q2a_ = UnitQuaternion(q2a.w, q2a.x, q2a.y, q2a.z, false)
        q2b_ = UnitQuaternion(q2b.w, q2b.x, q2b.y, q2b.z, false)
        velocity = A * (vrotate(ω1b,q2a_\q2b_) - ω1a) # in body1's frame
        force = -2 * Aᵀ * A * Diagonal(jr0.damper) * Aᵀ * velocity
        force = vrotate(force, q2b_ \ q2a_) # in body2's frame
        return force
    end
    ForwardDiff.jacobian(f, ω1b)
end

ω1a = rand(3)
q2a = UnitQuaternion(rand(4)...)
ω1b = rand(3)
q2b = UnitQuaternion(rand(4)...)
d1 = der1(ω1a, q2a, ω1b, q2b)
d3 = der3(ω1a, q2a, ω1b, q2b)
d4 = der4(ω1a, q2a, ω1b, q2b)
norm(d4 - d3)
norm(d4 - d1)

################################################################################
# Spring translation
################################################################################
jt0 = FixedOrientation(bodya0, bodyb0; spring = zeros(3), damper = zeros(3))[1][1]
jt0.spring = 1e1 .* rand(3)
jt0.damper = 1e1 .* rand(3)

# jt0 = Planar(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)
#
# jt0 = Prismatic(bodya0, bodyb0, rand(3); spring = zeros(3), damper = zeros(3))[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)
#
# jt0 = Fixed(bodya0, bodyb0)[1][1]
# jt0.spring = 1e1 .* rand(3)
# jt0.damper = 1e1 .* rand(3)

springforcea(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
springforceb(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
springforceb(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

Dspr1 = diagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Dspr2 = offdiagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Dspr3 = offdiagonal∂spring∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

fd_Dspr1 = fd_diagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dspr2 = fd_offdiagonal∂spring∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dspr3 = fd_offdiagonal∂spring∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

norm(Dspr1 - fd_Dspr1)
norm(Dspr2 - fd_Dspr2)
norm(Dspr3 - fd_Dspr3)
