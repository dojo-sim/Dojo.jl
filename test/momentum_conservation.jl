
vis = Visualizer()
open(vis)

################################################################################
# DICE
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
include("conservation_test.jl")

Δt0 = 0.01
g0 = 0.0
mech = getmechanism(:dice, Δt = Δt0, g = g0, contact = false)

v0 = [1,2,3.0]
ω0 = [1,1,1.0]
initialize!(mech, :dice, v = v0, ω = ω0)

storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:dice, Δt = Δt0, g = g0, contact = false)
    initialize!(mechanism, :dice, v = v0, ω = ω0)
    storage = simulate!(mechanism, t, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.0 + 0.2 * i for i = 1:20]
ms = getmomentum.(ts)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum" )
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum" )
@test all(norm.(ms, Inf) .< 1e-12)


################################################################################
# SINGLE PENDULUM
################################################################################
# single body
# initial angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
include("conservation_test.jl")

Δt0 = 0.01
g0 = 0.0
mech = getmechanism(:pendulum, Δt = Δt0, g = g0)

ϕ0 = 0.7
ω0 = 5.0
initialize!(mech, :pendulum, ϕ1 = ϕ0, ω1 = ω0)

storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:pendulum, Δt = Δt0, g = g0)
    initialize!(mechanism, :pendulum, ϕ1 = ϕ0, ω1 = ω0)
    storage = simulate!(mechanism, t, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.0 + 0.2 * i for i = 1:20]
ms = getmomentum.(ts)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum" )
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum" )
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-12)


################################################################################
# 2-lINK SNAKE SPHERICAL
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
include("conservation_test.jl")
Δt0 = 0.01
g0 = 0.0
Nlink0 = 2
spring0 = 0.0
damper0 = 0.0
mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Spherical, contact = false)

# ϕ0 = 0.7
# v0 = [-0.1,0.5,0.2]
# ω0 = [1,2,3.0]
# Δv0 = zeros(3)
# Δω0 = [3,-2,3.0] / Nlink0

ϕ0 = 0.7
v0 = 0.0*[-0.1,0.5,0.2]
ω0 = 1.0 * [2.0, 1.0, 3.0]
Δv0 = zeros(3)
Δω0 = 0.0 * [0,0,0.] / Nlink0
initialize!(mech, :snake, v = v0, ω = ω0)#, Δv = Δv0, Δω = Δω0)
# mech.bodies[3].J = Array(Diagonal([2.0, 2.9, 1.0]))
# mech.bodies[4].J = Array(Diagonal([1.0, 2.0, 3.0]))

storage = simulate!(mech, 25.0, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = :Spherical, contact = false)
    initialize!(mechanism, :snake, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)
    storage = simulate!(mechanism, t, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.0 + 0.2 * i for i = 1:20]
ms = getmomentum.(ts)
ms = [m .- ms[5] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum" )
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum" )
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11)
plot([q.w for q in storage.q[1]])
plot([q.w for q in storage.q[2]])
plot([w[1] for w in storage.ω[1]])
plot([w[2] for w in storage.ω[1]])
plot([w[3] for w in storage.ω[1]])
