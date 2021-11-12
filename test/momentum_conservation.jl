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

# Momentum computation
function getmomentum(model::Symbol, t::T, Δt::T, g::T, ϵ::T, controller!::Any; mech_kwargs::Dict = Dict(), init_kwargs::Dict = Dict()) where T
    mechanism = getmechanism(model, Δt = Δt, g = g; mech_kwargs...)
    initialize!(mechanism, model; init_kwargs...)
    storage = simulate!(mechanism, t, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ)
    # visualize(mechanism, storage, vis = vis)
    return DifferentiableContact.momentum(mechanism)
end


include(joinpath("..", "src", "optional_components", "energy.jl"))

# # visualizer
# vis = Visualizer()
# open(vis)

# Data
ϵ0 = 1e-14
Δt0 = 0.01
g0 = 0.0
ts = [0.5 + 0.2 * i for i = 1:5]
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

storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)

ms = getmomentum.(:dice, ts, Δt0, g0, ϵ0, controller!,
    init_kwargs = Dict(:v => v0, :ω => ω0))
ms = [m .- ms[1] for m in ms]
# plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
# plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@testset "Dice" begin @test all(norm.(ms, Inf) .< 1e-11) end


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

storage = simulate!(mech, 5.0, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)

ms = getmomentum.(:pendulum, ts, Δt0, g0, ϵ0, nocontrol!;
    init_kwargs = Dict(:ϕ1 => ϕ0, :ω1 => ω0))
ms = [m .- ms[1] for m in ms]
# plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
# plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@testset "Pendulum" begin @test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11) end


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
damper0 = 0.1
mech = getmechanism(:humanoid, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :humanoid)

storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false)
# visualize(mech, storage, vis = vis)

ms = getmomentum.(:humanoid, ts, Δt0, g0, ϵ0, controller!;
    mech_kwargs = Dict(:contact => false, :spring => spring0, :damper => damper0))
ms = [m .- ms[1] for m in ms]
plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@testset "Humanoid" begin @test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11) end


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
damper0 = 0.1
mech = getmechanism(:atlas, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :atlas)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)

ms = getmomentum.(:atlas, ts, Δt0, g0, ϵ0, controller!;
    mech_kwargs = Dict(:contact => false, :spring => spring0, :damper => damper0))
ms = [m .- ms[1] for m in ms]
# plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
# plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@testset "Atlas" begin @test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11) end


################################################################################
#  QUADRUPED
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 10.0
damper0 = 0.1
mech = getmechanism(:quadruped, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)

ms = getmomentum.(:quadruped, ts, Δt0, g0, ϵ0, controller!;
    mech_kwargs = Dict(:contact => false, :spring => spring0, :damper => damper0))
ms = [m .- ms[1] for m in ms]
# plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
# plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
@testset "Quadruped" begin @test all(norm.([m[4:6] for m in ms], Inf) .< 1e-11) end


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
damper0 = 1.0 * 4e0
mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Prismatic, contact = false)

v0 = 1.0 * [1, 2, 3] * Δt0
ω0 = 1.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.6*π))
initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, storage, vis = vis)

@testset "Snake" begin
    for jointtype in jointtypes
        ms = getmomentum.(:snake, ts, Δt0, g0, ϵ0, controller!;
            mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => jointtype),
            init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
        ms = [m .- ms[1] for m in ms]
        # plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
        # plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
        @test all(norm.(ms, Inf) .< 1e-11)
    end
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
damper0 = 1.0 * 4e-1
mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :PrismaticOrbital, contact = false)

v0 = 1.0 * [1, 2, 3] * Δt0
ω0 = 1.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.6*π))
initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# visualize(mech, downsample(storage, 1), vis = vis)

@testset "Twister" begin
    for jointtype in jointtypes
        ms = getmomentum.(:twister, ts, Δt0, g0, ϵ0, controller!;
            mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => jointtype),
            init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
        ms = [m .- ms[1] for m in ms]
        # plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
        # plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
        @test all(norm.(ms, Inf) .< 1e-11)
    end
end


#
#
# Nlink0 = 1
# spring0 = 0.0 * 4e0
# damper0 = 0.0 * 4e0
# mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
#     jointtype = :Prismatic, contact = false)
#
# v0 = 1.0 * [1, 2, 3] * Δt0
# ω0 = 100.0 * [1, 2, 3.0] * Δt0
# q10 = UnitQuaternion(RotX(0.6*π))
# initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
# storage = simulate!(mech, 5.0, controller!, record = true, solver = :mehrotra!, verbose = false, ϵ = ϵ0)
# # visualize(mech, storage, vis = vis)
#
# ts = [0.5 + 0.2 * i for i = 1:15]
# ms = getmomentum.(:snake, ts, Δt0, g0, ϵ0, controller!;
#     mech_kwargs = Dict(:Nlink => Nlink0, :contact => false, :spring => spring0, :damper => damper0, :jointtype => :Fixed),
#     init_kwargs = Dict(:q1 => q10, :v => v0, :ω => ω0))
# ms = [m .- ms[1] for m in ms]
# plot(ts, hcat(ms...)'[:,1:3], label = ["x" "y" "z"], title = "linear momentum")
# plot(ts, hcat(ms...)'[:,4:6], label = ["x" "y" "z"], title = "angular momentum")
# @test all(norm.(ms, Inf) .< 1e-11)
