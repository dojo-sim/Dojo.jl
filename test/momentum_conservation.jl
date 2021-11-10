
vis = Visualizer()
open(vis)

include(joinpath("..", "src", "optional_components", "energy.jl"))


function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = getcontroldim(joint)
        if nu <= 5
            if k ∈ (1:100)
                u = 0.5 * Δt0 * sones(nu)
            else
                u = szeros(nu)
            end
            setForce!(mechanism, joint, u)
        end
    end
    return
end

################################################################################
# DICE
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
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
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
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
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-12)


################################################################################
# 5-lINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
Nlink0 = 5
spring0 = 0.1
damper0 = 2e-1
mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Prismatic, contact = false)

v0 = 0.0 * [-0.1,0.5,0.2]
ω0 = 1.0 * [2.0, 1.0, 3.0]
Δv0 = zeros(3)
Δω0 = 0.0 * [1,2,-1.2] / Nlink0
initialize!(mech, :snake, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)


storage = simulate!(mech, 25.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T, jointtype::Symbol) where T
    mechanism = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = jointtype, contact = false)
    initialize!(mechanism, :snake, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

for jointtype in (:Revolute, :Orbital, :Spherical, :Fixed)#, :Prismatic, :Planar, :FixedOrientation)
# for jointtype in (:Prismatic, :Planar)#, :FixedOrientation, :Fixed)
    ts = [1.0 + 0.2 * i for i = 1:10]
    ms = getmomentum.(ts, jointtype)
    ms = [m .- ms[1] for m in ms]
    plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
    plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
    @test all(norm.(ms, Inf) .< 1e-11)
end


ts = [1.0 + 0.2 * i for i = 1:10]
ms = getmomentum.(ts, :Prismatic)
# ms = getmomentum.(ts, :Prismatic)
# ms = getmomentum.(ts, :FixedOrientation)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.(ms, Inf) .< 1e-11)

################################################################################
#  HUMANOID
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
spring0 = 0.1
damper0 = 0.1
mech = getmechanism(:humanoid, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :humanoid)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:humanoid, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
    initialize!(mechanism, :humanoid)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.5 + 0.2 * i for i = 1:5]
ms = getmomentum.(ts)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11)


################################################################################
#  ATLAS
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
spring0 = 0.1
damper0 = 0.1
mech = getmechanism(:atlas, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :atlas)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:atlas, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
    initialize!(mechanism, :atlas)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.5 + 0.2 * i for i = 1:5]
ms = getmomentum.(ts)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11)


################################################################################
#  QUADRUPED
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
spring0 = 0.1
damper0 = 0.1
mech = getmechanism(:quadruped, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T) where T
    mechanism = getmechanism(:atlas, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
    initialize!(mechanism, :atlas)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.5 + 0.2 * i for i = 1:5]
ms = getmomentum.(ts)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11)






















################################################################################
# 5-lINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
Nlink0 = 2
spring0 = 0.0 * 3e0
damper0 = 0.0 * 2e-1
mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Planar, contact = false)

eqc2 = collect(mech.eqconstraints)[2]
tra2 = eqc2.constraints[1]
tra2.spring
tra2.damper
rot2 = eqc2.constraints[2]
rot2.spring
rot2.damper

eqc2.isspring
eqc2.isdamper

v0 = 0.0 * [-0.1,0.5,0.2]
ω0 = 1.0 * [2.0, 1.0, 3.0]
Δv0 = zeros(3)
Δω0 = 0.0 * [1,2,-1.2] / Nlink0
initialize!(mech, :snake, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = getcontroldim(joint)
        @show nu
        if nu <= 5
            if k ∈ (1:100)
                u = -1.0 * Δt0 * sones(nu)
            else
                u = szeros(nu)
            end
            setForce!(mechanism, joint, u)
        end
    end
    return
end


storage = simulate!(mech, 3.3, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, storage, vis = vis)

function getmomentum(t::T, jointtype::Symbol) where T
    mechanism = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = jointtype, contact = false)
    initialize!(mechanism, :snake, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

# for jointtype in (:Revolute, :Orbital, :Spherical, :Fixed)#, :Prismatic, :Planar, :FixedOrientation)

ts = [1.0 + 0.2 * i for i = 1:10]
ms = getmomentum.(ts, :Planar)
# ms = getmomentum.(ts, :Prismatic)
# ms = getmomentum.(ts, :Prismatic)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.(ms, Inf) .< 1e-11)














vis = Visualizer()
open(vis)


################################################################################
# 5-lINK TWISTER
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Δt0 = 0.01
g0 = 0.0
Nlink0 = 2
spring0 = 1.0 * 3e0
damper0 = 1.0 * 3e0
mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Prismatic, contact = false)

v0 = 1.0 * [25π, 0, 0.0] * Δt0
ω0 = 1.0 * [0, 0, 0.0] * Δt0
Δv0 = 1.0 * [0, 0, 0.0] * Δt0
Δω0 = 1.0 * [0, 50*2π, 0.0] * Δt0
q10x = UnitQuaternion(RotY(π/2))
q10y = UnitQuaternion(RotX(π/2))
q10 = UnitQuaternion(RotX(π))

initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)

function controller!(mechanism, k)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = getcontroldim(joint)
        if nu <= 5
            if k ∈ (1:100)
                u = 0.0 * Δt0 * sones(nu)
            else
                u = szeros(nu)
            end
            setForce!(mechanism, joint, u)
        end
    end
    return
end

storage = simulate!(mech, 550, controller!, record = true, solver = :mehrotra!, verbose = false)
visualize(mech, downsample(storage, 10), vis = vis)
plot([ω[2] for ω in storage.ω[1]])

function getmomentum(t::T, jointtype::Symbol) where T
    mechanism = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
        jointtype = jointtype, contact = false)
    initialize!(mechanism, :twister, q1 = q10, v = v0, ω = ω0, Δv = Δv0, Δω = Δω0)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false)
    return momentum(mechanism)
end

ts = [1.0 + 5.2 * i for i = 1:10]
ms = getmomentum.(ts, :Prismatic)
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@test all(norm.(ms, Inf) .< 1e-11)

eqc2 = collect(mech.eqconstraints)[2]
tra2 = eqc2.constraints[1]
tra2.spring
tra2.damper
rot2 = eqc2.constraints[2]
rot2.spring
rot2.damper

eqc2.isspring
eqc2.isdamper

eqcs2.axis
zerodimstaticadjoint(constraintmat(tra2)) * constraintmat(tra2)
constraintmat(rot2)