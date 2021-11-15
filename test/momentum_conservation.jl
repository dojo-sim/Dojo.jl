# visualizer
vis = Visualizer()
open(vis)

const DifferentiableContact = Main

# Controller
function controller!(mechanism, k; U = 0.5, Δt = 0.01)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = getcontroldim(joint)
        u = (nu <= 5 && k ∈ (1:100)) * U * Δt * sones(nu)
        setForce!(mechanism, joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U = 0.0)


# Data
ϵ0 = 1e-14
Δt0 = 0.01
g0 = 0.0
jointtypes = [
    :Fixed,
    :Prismatic,
    :Planar,
    :FixedOrientation,
    :Revolute,
    :Cylindrical,
    :PlanarAxis,
    :FreeRevolute,
    :Orbital,
    :PrismaticOrbital,
    :PlanarOrbital,
    :FreeOrbital,
    :Spherical,
    :CylindricalFree,
    :PlanarFree
    ]

################################################################################
# DICE
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
mech = getmechanism(:dice, Δt = Δt0, g = g0, contact = false)
v0 = [1,2,3.0]
ω0 = [1,1,1.0]
initialize!(mech, :dice, v = v0, ω = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
ke0 = kineticEnergy(mech, storage)[5:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
@test norm(ke0 .- ke0[1], Inf) < 1e-11

################################################################################
# SINGLE PENDULUM
################################################################################
# single body
# initial angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
mech = getmechanism(:pendulum, Δt = Δt0, g = g0)
ϕ0 = 0.7
ω0 = 5.0
initialize!(mech, :pendulum, ϕ1 = ϕ0, ω1 = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[10:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
ke0 = kineticEnergy(mech, storage)[10:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
# @test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
@test norm(ke0 .- ke0[1], Inf) < 1e-11



################################################################################
#  HUMANOID
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 10.0
damper0 = 1.0
mech = getmechanism(:humanoid, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :humanoid)
bodies = collect(mech.bodies)
setVelocity!.(bodies, ω = 1e-0rand(3))


storage = simulate!(mech, 10.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, downsample(storage, 10), vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# ke0 = kineticEnergy(mech, storage)[1:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11


################################################################################
#  ATLAS
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 10.0
damper0 = 1.0
mech = getmechanism(:atlas, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :atlas)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# ke0 = kineticEnergy(mech, storage)[1:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11



################################################################################
#  QUADRUPED
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 0.3
damper0 = 0.1
mech = getmechanism(:quadruped, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, controller_quadruped!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# ke0 = kineticEnergy(mech, storage)[1:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11




################################################################################
# 5-lINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Nlink0 = 5
spring0 = 1.0 * 4e0
damper0 = 1.0 * 2e+1


mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Revolute, contact = false, r = 0.05)

v0 = 100.0 * [1, 2, 3] * Δt0
ω0 = 100.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 1.50, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)


m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
ke0 = kineticEnergy(mech, storage)[5:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11

for jointtype in jointtypes
    @show jointtype
    mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = jointtype, contact = false, r = 0.05)

    v0 = 100.0 * [1, 2, 3] * Δt0
    ω0 = 100.0 * [1, 2, 3.0] * Δt0
    q10 = UnitQuaternion(RotX(0.5*π))
    initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
    storage = simulate!(mech, 1.50, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
    visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[5:end]
    mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
    mang0 = [Vector(m-m0[1])[4:6] for m in m0]
    plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
    plt = plot()
    plot!([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
    display(plt)
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

################################################################################
# 5-lINK TWISTER
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Nlink0 = 5
spring0 = 1.0 * 4e0
damper0 = 1.0 * 2e+1

mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Revolute, contact = false, r = 0.05)

v0 = 100.0 * [1, 2, 3] * Δt0
ω0 = 100.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 1.50, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
visualize(mech, storage, vis = vis)


m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
ke0 = kineticEnergy(mech, storage)[5:end]
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
@test all(norm.(mlin0, Inf) .< 1e-11)
@test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11

for jointtype in jointtypes
    @show jointtype
    mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = jointtype, contact = false, r = 0.05)

    v0 = 100.0 * [1, 2, 3] * Δt0
    ω0 = 100.0 * [1, 2, 3.0] * Δt0
    q10 = UnitQuaternion(RotX(0.5*π))
    initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
    storage = simulate!(mech, 1.50, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
    visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[5:end]
    mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
    mang0 = [Vector(m-m0[1])[4:6] for m in m0]
    plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
    plt = plot()
    plot!([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
    display(plt)
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end




































#
#
#
#
# kineticEnergy(mech)
# body1 = collect(mech.bodies)[1]
# body2 = collect(mech.bodies)[2]
# 0.5 * body1.m * v0'*v0 + 0.5 * body1.m * v0'*v0
# Δt0 = 1e-3
# ts = [0.1 + 0.1 * i for i = 1:15]
# g0 = 0.0
# ϵ0 = 1e-14
# ms = getenergy.(:snake, ts, Δt0, g0, ϵ0, nocontrol!; complete = true,
#     mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => :Fixed, :r => 0.005),
#     init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
#
# plot(ts, [m[1] for m in ms], label = "energy", title = "mechanical energy", linewidth = 5)
# plot!(ts, [m[2] for m in ms], label = "energy", title = "potential energy", linewidth = 5)
# plot!(ts, [m[3] for m in ms], label = "energy", title = "kinetic energy", linewidth = 5)
# plot(ts, [m[1]-ms[1][1] for m in ms], label = "energy", title = "mechanical energy", linewidth = 5)
# plot(ts, [m[2]-ms[1][2] for m in ms], label = "energy", title = "potential energy", linewidth = 5)
# plot(ts, [m[3]-ms[1][3] for m in ms], label = "energy", title = "kinetic energy", linewidth = 5)
# plot(ts, hcat([m[4:6] .- ms[1][4:6] for m in ms]...)', label = ["x" "y" "z"], title = "linear momentum", linewidth = 5)
# plot(ts, hcat([m[7:9] .- ms[1][7:9] for m in ms]...)', label = ["x" "y" "z"], title = "angular momentum", linewidth = 5)
# body1.state.τk
# body1.state.Fk
#
#
# body1 = collect(mech.bodies)[1]
# body1.state.ωsol[2] - body1.state.ωsol[1]
#
#
# ################################################################################
# # 5-lINK SNAKE
# ################################################################################
# # multiple bodies
# # initial linear and angular velocity
# # no gravity
# # with spring and damper
# # with control
# ################################################################################
# Nlink0 = 1
# Δt0 = 0.01
# spring0 = 0.0 * 4e0
# damper0 = 0.0 * 4e0
# mech1 = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :Fixed, contact = false, h = 2.0)
#
# # v0 = 1.0 * [1, 2, 3] * Δt0
# # ω0 = 100.0 * [1, 2, 3.0] * Δt0
#
# # v0 = 400.0 * [0, 0, -0/2] * Δt0
# # ω0 = 400.0 * [1, 0, 0.0] * Δt0
#
# # v0 = 400.0 * [0/2, 0, 0.0] * Δt0
# # ω0 = 400.0 * [0, 1.0, 0.0] * Δt0
#
# v0 = 400.0 * [0/2, 0/2, -0/2] * Δt0
# ω0 = 400.0 * [1, 1.0, 1.0] * Δt0
#
# q10 = UnitQuaternion(RotX(0.5*π))
# initialize!(mech1, :snake, q1 = q10, v = v0, ω = ω0)
# storage1 = simulate!(mech1, 6.0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech1, storage1, vis = vis)
#
# kineticEnergy(mech, linear = true, angular = false)
# kineticEnergy(mech1, linear = true, angular = false)
#
#
#
#
# ts = [0.1 + 0.1 * i for i = 1:20]
# ms = getenergy.(:snake, ts, Δt0, g0, ϵ0, controller!; complete = true,
#     mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => :Fixed, :h => 2.0),
#     init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
#
# plot(ts, [m[1] for m in ms], label = "energy", title = "mechanical energy", linewidth = 5)
# plot!(ts, [m[2] for m in ms], label = "energy", title = "potential energy", linewidth = 5)
# plot!(ts, [m[3] for m in ms], label = "energy", title = "kinetic energy", linewidth = 5)
# plot(ts, [m[1]-ms[1][1] for m in ms], label = "energy", title = "mechanical energy", linewidth = 5)
# plot(ts, [m[2]-ms[1][2] for m in ms], label = "energy", title = "potential energy", linewidth = 5)
# plot(ts, [m[3]-ms[1][3] for m in ms], label = "energy", title = "kinetic energy", linewidth = 5)
# plot(ts, hcat([m[4:6] .- ms[1][4:6] for m in ms]...)', label = ["x" "y" "z"], title = "linear momentum", linewidth = 5)
# plot(ts, hcat([m[7:9] .- ms[1][7:9] for m in ms]...)', label = ["x" "y" "z"], title = "angular momentum", linewidth = 5)
#
#
#
#
#
#
#
#
# tsim0 = 0.005
# Δt0 = 1e-2
# α = 400
# Nlink0 = 2
#
# mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :Fixed, contact = false, r = 0.005)
#
# # v0 = 400.0 * [1/2, 0, -1/2] * Δt0
# # ω0 = 400.0 * [1, 1.0, 0.0] * Δt0
# v0 = α * [1/2, 0, -1/2] * Δt0
# ω0 = α * [1, 1.0, 1.0] * Δt0
#
# initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
# storage = simulate!(mech, tsim0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)
#
# storage.q[1][1]
# storage.q[2][1]
#
# storage.q[1][2]
# storage.q[2][2]
#
# storage.v[1][1]
# storage.v[2][1]
#
# storage.ω[1][1]
# storage.ω[2][1]
#
#
#
# Nlink0 = 1
# mech1 = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :Fixed, contact = false, h = 2.0, r = 0.005)
#
# # v0 = 400.0 * [0/2, 0, -0/2] * Δt0
# # ω0 = 400.0 * [1, 1.0, 0.0] * Δt0
# v0 = α * [0/2, 0, -0/2] * Δt0
# ω0 = α * [1, 1.0, 1.0] * Δt0
#
# initialize!(mech1, :snake, q1 = q10, v = v0, ω = ω0)
# storage1 = simulate!(mech1, tsim0, nocontrol!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech1, storage1, vis = vis)
#
# kineticEnergy(mech, linear = true, angular = true) - kineticEnergy(mech1, linear = true, angular = true)
# kineticEnergy(mech, linear = true, angular = true)
# kineticEnergy(mech1, linear = true, angular = true)
#
# body11 = collect(mech1.bodies)[1]
# body21 = collect(mech.bodies)[1]
# body22 = collect(mech.bodies)[2]
#
# q11 = posargsk(body11.state)[2]
# plin11 = momentum_body(mech1, body11)[1:3]
# # pang11 = rotation_matrix(inv(q11)) * momentum_body(mech1, body11)[4:6]
# pang11 = rotation_matrix(inv(q21)) * momentum_body(mech1, body11)[4:6]
#
# momentum_body(mech1, body11)[4:6] ./ momentum_body(mech, body21)[4:6]
# momentum_body(mech, body22)[4:6]
#
# rotation_matrix(q11) - rotation_matrix(q21)
# rotation_matrix(q11) - rotation_matrix(q22)
# rotation_matrix(q21) - rotation_matrix(q22)
#
# [q11.w, q11.x, q11.y, q11.z] - [q21.w, q21.x, q21.y, q21.z]
# [q22.w, q22.x, q22.y, q22.z] - [q21.w, q21.x, q21.y, q21.z]
#
# ET1 = 0.5 * pang11' * (body11.J \ pang11)
#
# q21 = posargsk(body21.state)[2]
# q22 = posargsk(body22.state)[2]
# plin21 = momentum_body(mech, body21)[1:3]
# plin22 = momentum_body(mech, body22)[1:3]
# v21 = plin21/body21.m
# v22 = plin22/body22.m
# ω11 = body11.J \ pang11
# norm(skew(ω11) * [0,0,0.5] * body21.m)
# norm(plin21)
#
#
# pang21 = rotation_matrix(inv(q21)) * momentum_body(mech, body21)[4:6]
# pang22 = rotation_matrix(inv(q22)) * momentum_body(mech, body22)[4:6]
# ET2 = 0.5 * pang21' * (body21.J \ pang21) +
#     0.5 * pang22' * (body22.J \ pang22) +
#     0.5 * plin21' * (body21.m \ plin21) +
#     0.5 * plin22' * (body22.m \ plin22)
# 0.5 * pang21' * (body21.J \ pang21)
# 0.5 * pang22' * (body22.J \ pang22)
# 0.5 * plin21' * (body21.m \ plin21)
# 0.5 * plin22' * (body22.m \ plin22)
#
#
# ET1/(0.5 * pang21' * (body21.J \ pang21))
#
# body11.m
# body21.m
# body22.m
# body11.J
# body21.J
# body22.J
# r = 0.005
# h1 = 2.0
# h2 = 1.0
# norm(1/12 * body11.m * Diagonal([4r^2 + h1^2, 9r^2 + h1^2, 4r^2 + 9r^2]) - body11.J)
# norm(1/12 * body21.m * Diagonal([4r^2 + h2^2, 9r^2 + h2^2, 4r^2 + 9r^2]) - body21.J)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#







#
#
# ################################################################################
# # 5-lINK SNAKE
# ################################################################################
# # multiple bodies
# # initial linear and angular velocity
# # no gravity
# # with spring and damper
# # with control
# ################################################################################
# Nlink0 = 5
# spring0 = 1.0 * 4e0
# damper0 = 1.0 * 4e0
#
# v0 = 1.0 * [1, 2, 3] * Δt0
# ω0 = 1.0 * [1, 2, 3.0] * Δt0
# q10 = UnitQuaternion(RotX(0.6*π))
#
# mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :Prismatic, contact = false)
# initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
# storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)
#
# m0 = momentum(mech, storage)[1:end]
# mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
# mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# ke0 = kineticEnergy(mech, storage)[1:end]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
# @test all(norm.(mlin0, Inf) .< 1e-11)
# @test all(norm.(mang0, Inf) .< 1e-11)
# @test norm(ke0 .- ke0[1], Inf) < 1e-11
#
# @testset "Snake" begin
#     for jointtype in jointtypes
#         mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#             jointtype = jointtype, contact = false)
#         initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
#         storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
#         # visualize(mech, storage, vis = vis)
#
#         m0 = momentum(mech, storage)[1:end]
#         mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
#         mang0 = [Vector(m-m0[1])[4:6] for m in m0]
#         ke0 = kineticEnergy(mech, storage)[1:end]
#         # plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
#         # plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
#         # plot([(i-1)*Δt0 for i in 1:length(m0)], ke0 .- ke0[1])
#         @test all(norm.(mlin0, Inf) .< 1e-11)
#         @test all(norm.(mang0, Inf) .< 1e-11)
#         # @test norm(ke0 .- ke0[1], Inf) < 1e-11
#     end
# end
#
# ################################################################################
# # 5-lINK TWISTER
# ################################################################################
# # multiple bodies
# # initial linear and angular velocity
# # no gravity
# # with spring and damper
# # with control
# ################################################################################
# Nlink0 = 5
# spring0 = 1.0 * 4e0
# damper0 = 1.0 * 4e-1
# mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :PrismaticOrbital, contact = false)
#
# v0 = 1.0 * [1, 2, 3] * Δt0
# ω0 = 1.0 * [1, 2, 3.0] * Δt0
# q10 = UnitQuaternion(RotX(0.6*π))
# initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
# storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# # visualize(mech, downsample(storage, 1), vis = vis)
#
# @testset "Twister" begin
#     for jointtype in jointtypes
#         ms = getenergy.(:twister, ts, Δt0, g0, ϵ0;
#             mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => jointtype),
#             init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
#         ms = ms .- ms[1]
#         plot(ts, ms, label = "energy", title = "mechanical energy")
#         @test all(norm.(ms, Inf) .< 1e-11)
#     end
# end
