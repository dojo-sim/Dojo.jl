
# visualizer
# vis = Visualizer()
# open(vis)

# Controller
function controller!(mechanism, k; U = 0.5, Δt = 0.01)
    for (i,joint) in enumerate(mechanism.eqconstraints)
        nu = controldim(joint)
        u = (nu <= 5 && k ∈ (1:100)) * U * Δt * sones(nu)
        setForce!(mechanism, joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U = 0.0)


# Data
ϵ0 = 1e-12
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
# Box
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
mech = getmechanism(:box, Δt = Δt0, g = g0, contact = false)
v0 = [1,2,3.0]
ω0 = [1,1,1.0]
initialize!(mech, :box, v = v0, ω = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Box" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

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
ω0 = 3.0
initialize!(mech, :pendulum, ϕ1 = ϕ0, ω1 = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[10:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
# @test all(norm.(mlin0, Inf) .< 1e-11)
@testset "Momentum: Pendulum" begin
    @test all(norm.(mang0, Inf) .< 1e-11)
end

################################################################################
#  HUMANOID
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 0.0
damper0 = 1.0
mech = getmechanism(:humanoid, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :humanoid)
bodies = collect(mech.bodies)
setVelocity!.(bodies, ω = 1e-1 * rand(3))


storage = simulate!(mech, 10.0, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, downsample(storage, 1), vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Humanoid" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

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
storage = simulate!(mech, 5.0, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Atlas" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

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
damper0 = 1.0
mech = getmechanism(:quadruped, Δt = Δt0, g = g0, spring = spring0, damper = damper0, contact = false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 1.0, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Quadruped" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

################################################################################
# 5-lINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
Nlink0 = 2
spring0 = 1.0 * 4e0
damper0 = 1.0 * 2e+1

mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :Revolute, contact = false, r = 0.05);

v0 = 100.0 * [1, 2, 3] * Δt0
ω0 = 10.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 1.50, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Snake" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

@testset "Momentum: Snake" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = getmechanism(:snake, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 100.0 * [1, 2, 3] * Δt0
        ω0 = 10.0 * [1, 2, 3.0] * Δt0
        q10 = UnitQuaternion(RotX(0.5*π))
        initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 1.50, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
        mang0 = [Vector(m-m0[1])[4:6] for m in m0]
        # plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
        # plt = plot()
        # plot!([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
        # display(plt)
        @test all(norm.(mlin0, Inf) .< 1e-11)
        @test all(norm.(mang0, Inf) .< 1e-11)
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
damper0 = 1.0 * 2e+1

mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
    jointtype = :FixedOrientation, contact = false, r = 0.05)

v0 = 10.0 * [1, 2, 3] * Δt0
ω0 = 1.0 * [1, 2, 3.0] * Δt0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 1.25, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0 = [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Twister" begin
    @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1e-11)
end

@testset "Momentum: Twister" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = getmechanism(:twister, Δt = Δt0, g = g0, Nlink = Nlink0, spring = spring0, damper = damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 10.0 * [1, 2, 3] * Δt0
        ω0 = 1.0 * [1, 2, 3.0] * Δt0
        q10 = UnitQuaternion(RotX(0.5*π))
        initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 1.50, controller!, record = true, verbose = false, opts=InteriorPointOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
        mang0 = [Vector(m-m0[1])[4:6] for m in m0]
        # plot([(i-1)*Δt0 for i in 1:length(m0)], hcat(mlin0...)')
        # plt = plot()
        # plot!([(i-1)*Δt0 for i in 1:length(m0)], hcat(mang0...)')
        # display(plt)
        @test all(norm.(mlin0, Inf) .< 1e-11)
        @test all(norm.(mang0, Inf) .< 1e-11)
    end
end
