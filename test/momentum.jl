
# visualizer
# vis = Visualizer()
# open(vis)

# const Dojo = Main

# Controller
function controller!(mechanism, k; U = 0.5, timestep = 0.01)
    for (i,joint) in enumerate(mechanism.joints)
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:100)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U = 0.0)


# Data
ϵ0 = 1e-12
timestep0 = 0.01
gravity0= 0.0
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
mech = get_mechanism(:box, timestep=timestep0, gravity=gravity0, contact = false)
v0 = [1,2,3.0]
ω0 = [10,10,10.0]
initialize!(mech, :box, v = v0, ω = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
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
mech = get_mechanism(:pendulum, timestep=timestep0, gravity=gravity0)
ϕ0 = 0.7
ω0 = 5.0
initialize!(mech, :pendulum, ϕ1 = ϕ0, ω1 = ω0)

storage = simulate!(mech, 5.0, nocontrol!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[10:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
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
spring0 = 1.0
damper0 = 1.0
mech = get_mechanism(:humanoid, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0, contact = false)
initialize!(mech, :humanoid)
bodies = collect(mech.bodies)
set_maximal_velocities!.(bodies, ω = 1e-0rand(3))

storage = simulate!(mech, 10.0, controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Humanoid" begin
    @test all(norm.(mlin0, Inf) .< 1e-8)
    @test all(norm.(mang0, Inf) .< 1e-8)
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
mech = get_mechanism(:atlas, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0, contact = false)
initialize!(mech, :atlas)
storage = simulate!(mech, 5.0, controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Atlas" begin
    @test all(norm.(mlin0, Inf) .< 1e-8)
    @test all(norm.(mang0, Inf) .< 1e-8)
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
damper0 = 0.1
mech = get_mechanism(:quadruped, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0, contact = false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[1:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Quadruped" begin
    @test all(norm.(mlin0, Inf) .< 1e-8)
    @test all(norm.(mang0, Inf) .< 1e-8)
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
Nb0 = 5
spring0 = 0.0 * 4e0
damper0 = 0.0 * 2e+1

mech = get_mechanism(:snake, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
    jointtype = :Revolute, contact = false, r = 0.05);

v0 = 100.0 * [1, 2, 3] * timestep0
ω0 = 100.0 * [1, 2, 3.0] * timestep0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 1.50, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Snake" begin
    @test all(norm.(mlin0, Inf) .< 1e-8)
    @test all(norm.(mang0, Inf) .< 1e-8)
end

@testset "Momentum: Snake" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = get_mechanism(:snake, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3.0] * timestep0
        q10 = UnitQuaternion(RotX(0.5*π))
        initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 1.50, controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        # @show mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
        # @show mang0= [Vector(m-m0[1])[4:6] for m in m0]
        # plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
        # plt = plot()
        # plot!([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
        # display(plt)
        @test all(norm.(mlin0, Inf) .< 1e-8)
        @test all(norm.(mang0, Inf) .< 1e-8)
    end
end

################################################################################
# 5-lINK TWISTER
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
################################################################################

Nb0 = 5
spring0 = 1.0 * 4e0
damper0 = 1.0 * 2e+1

mech = get_mechanism(:twister, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
    jointtype = :FixedOrientation, contact = false, r = 0.05);

v0 = 100.0 * [1, 2, 3] * timestep0
ω0 = 100.0 * [1, 2, 3.0] * timestep0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 2.5, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

m0 = momentum(mech, storage)[5:end]
mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
mang0= [Vector(m-m0[1])[4:6] for m in m0]
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mlin0...)')
# plot([(i-1)*timestep0 for i in 1:length(m0)], hcat(mang0...)')
@testset "Momentum: Twister" begin
    @test all(norm.(mlin0, Inf) .< 1e-8)
    @test all(norm.(mang0, Inf) .< 1e-8)
end

@testset "Momentum: Twister" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = get_mechanism(:twister, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3.0] * timestep0
        q10 = UnitQuaternion(RotX(0.5*π))
        initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 1.50, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        mlin0 = [Vector(m-m0[1])[1:3] for m in m0]
        mang0= [Vector(m-m0[1])[4:6] for m in m0]

        @test all(norm.(mlin0, Inf) .< 1e-8)
        @test all(norm.(mang0, Inf) .< 1e-8)
    end
end