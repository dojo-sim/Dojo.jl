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
initialize!(mech, :npendulum, ϕ1 = 0.7)

j1 = mech.eqconstraints[1]
jt1 = j1.constraints[1]
jr1 = j1.constraints[2]
j1.isdamper = true
j1.isspring = true

jr1.spring = 1e1 .* sones(3)# 1e4
jr1.damper = 1e1 .* sones(3)# 1e4
mech.eqconstraints[1].isdamper
mech.eqconstraints[1].constraints[2].damper

storage = simulate!(mech, 0.1, record = true, solver = :mehrotra!)
# plot(hcat(Vector.(storage.x[1])...)')
# plot(hcat([[q.w, q.x, q.y, q.z] for q in storage.q[1]]...)')
# plot(hcat(Vector.(storage.v[1])...)')
# plot(hcat(Vector.(storage.ω[1])...)')

# visualize(mech, storage, vis = vis)

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

Dtra1 = diagonal∂damper∂ʳvel(jt0)
Dtra2 = offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
Dtra3 = offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

fd_Dtra1 = fd_diagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dtra2 = fd_offdiagonal∂damper∂ʳvel(jt0, x2a0, q2a0, x2b0, q2b0, x1a0, v1a0, q1a0, ω1a0, x1b0, v1b0, q1b0, ω1b0, Δt0)
fd_Dtra3 = fd_offdiagonal∂damper∂ʳvel(jt0, x2b0, q2b0, x1b0, v1b0, q1b0, ω1b0, Δt0)

norm(Dtra1 - fd_Dtra1)
norm(Dtra2 + fd_Dtra2)
norm(Dtra3 + fd_Dtra3)



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

Drot1 = diagonal∂damper∂ʳvel(jr0)
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
