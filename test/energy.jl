# # visualizer
# vis = Visualizer()
# open(vis)

# const Dojo = Main

# Data
ϵ0 = 1e-12
timestep0 = 1e-2
start0 = Int(floor(1/timestep0)) + 1
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

# Controller
function controller!(mechanism, k; U=0.5, timestep=timestep0)
    N = Int(floor(1 / timestep))
    for (i,joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U = 0.0)

################################################################################
# DICE
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
gravity0 = -10.0
mech = get_mechanism(:box, timestep=timestep0, gravity=gravity0, contact=false)
v0 = [1,2,3.0]
ω0 = [1,1,1.0]
initialize!(mech, :box, v=v0, ω=ω0)

storage = simulate!(mech, 5.0, nocontrol!, record=true, verbose=false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
# With spring the error is or the order of timestep, but there is no drift
@testset "Energy: Dice" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-11
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)


################################################################################
# SINGLE PENDULUM
################################################################################
# single body
# initial angular velocity
# gravity
# no damper
# no control
################################################################################
gravity0 = -10.0
spring0 = 1.0
damper0 = 0.0
mech = get_mechanism(:pendulum, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0)
ϕ0 = 0.9π
ω0 = 2π

initialize!(mech, :pendulum, ϕ1=ϕ0, ω1=ω0)
storage = simulate!(mech, 5.0, controller!, record=true, verbose=false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
# With spring the error is or the order of timestep, but there is no drift
@testset "Energy: Pendulum" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

################################################################################
# SLIDER
################################################################################
# single body
# no initial linear velocity
# no gravity
# no damper
# no control
################################################################################
gravity0 = -0.0
spring0 = 10.0
damper0 = 0.0
mech = get_mechanism(:slider, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0)
z0 = 0.5
initialize!(mech, :slider, z1=z0)

# Analytical
body1 = collect(mech.bodies)[1]
zmax = z0
vmax = z0 * sqrt(spring0 / body1.mass)
pe_max = 0.5 * spring0 * zmax^2
ke_max = 0.5 * body1.mass * vmax^2

storage = simulate!(mech, 5.0,  nocontrol!, record=true, verbose=false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

@testset "Energy: Slider" begin
    # Test maximum amplitude and velocity
    @test norm(maximum([x[3] for x in storage.x[1]]) - zmax + 0.5) < 1e-4
    @test norm(maximum([vl[3] for vl in storage.vl[1]]) - vmax) < 1e-4

    # Test mechanical energy conservation
    # With spring the error is or the order of timestep, but there is no drift
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
    # @show norm((me0 .- me0[1]) ./ mean(me0), Inf)
end

################################################################################
# SLIDER
################################################################################
# single body
# no initial linear velocity
# gravity
# no damper
# no control
################################################################################
gravity0 = -9.81
spring0 = 0.0
damper0 = 0.0
mech = get_mechanism(:slider, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0)
z0 = 0.5
initialize!(mech, :slider, z1 = z0)

storage = simulate!(mech, 1.5,  nocontrol!, record=true, verbose=false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
# For gravity the conservation is perfect
@testset "Energy: Slider" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
end
norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-9


################################################################################
# SLIDER
################################################################################
# single body
# no initial linear velocity
# gravity
# no damper
# no control
################################################################################
gravity0 = -9.81
spring0 = 10.0
damper0 = 0.0
mech = get_mechanism(:slider, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0)
z0 = 0.5
initialize!(mech, :slider, z1 = z0)

storage = simulate!(mech, 5.0,  nocontrol!, record=true, verbose=false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[start0])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[start0])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[start0])

# Test mechanical energy conservation
@testset "Energy: Slider" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)


################################################################################
#  HUMANOID
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
function humanoid_controller!(mechanism, k; U=0.05, timestep=timestep0)
    N = Int(floor(1/timestep))
    for (i,joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

gravity0 = 0.0
spring0 = 1.0
damper0 = 0.0
mech = get_mechanism(:humanoid, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0, contact = false)
initialize!(mech, :humanoid)
bodies = collect(mech.bodies)
for body in mech.bodies
    set_velocity!(body, ω = 0.5*rand(3))
    # set_velocity!(body, v = 1.0*rand(3))
end

storage = simulate!(mech, 3.0, humanoid_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, downsample(storage, 1), vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Humanoid" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 2e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

################################################################################
#  ATLAS
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
spring0 = 1.0
damper0 = 0.0
mech = get_mechanism(:atlas, timestep=timestep0, gravity=gravity0, spring=spring0, damper=damper0, contact = false)
initialize!(mech, :atlas)
bodies = collect(mech.bodies)
set_velocity!.(bodies, ω = 1.0*rand(3))

storage = simulate!(mech, 5.0, humanoid_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Atlas" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 3e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)



################################################################################
#  QUADRUPED
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
function quadruped_controller!(mechanism, k; U = 0.01, timestep=timestep0)
    N = Int(floor(1/timestep))
    for (i,joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

spring0 = 0.1
damper0 = 0.0
mech = get_mechanism(:quadruped, timestep=timestep0, gravity=gravity0, spring=spring0,
    damper=damper0, contact = false, limits=false)
initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, quadruped_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Quadruped" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)


################################################################################
# 5-lINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
function snake_controller!(mechanism, k; U = 0.05, timestep=timestep0)
    N = Int(floor(1/timestep))
    for (i,joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

Nb0 = 5
spring0 = 0.01
damper0 = 0.0
mech = get_mechanism(:snake, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
    jointtype = :Revolute, contact = false, r = 0.05);

v0 = 10.0 * [1, 2, 3] * timestep0
ω0 = 10.0 * [1, 2, 3] * timestep0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 3.0, snake_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Snake" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

@testset "Energy: Snake" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = get_mechanism(:snake, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3] * timestep0
        q10 = UnitQuaternion(RotX(0.5*π))

        initialize!(mech, :snake, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 3.0, snake_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        ke0 = kinetic_energy(mech, storage)[start0:end]
        pe0 = potential_energy(mech, storage)[start0:end]
        me0 = mechanical_energy(mech, storage)[start0:end]

        # plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
        # plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
        # plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

        # Test mechanical energy conservation
        @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
        norm((me0 .- me0[1]) ./ mean(me0), Inf)
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
Nb0 = 5
spring0 = 0.01
damper0 = 0.0
mech = get_mechanism(:twister, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
    jointtype = :Revolute, contact = false, r = 0.05);

v0 = 10.0 * [1, 2, 3] * timestep0
ω0 = 10.0 * [1, 2, 3] * timestep0
q10 = UnitQuaternion(RotX(0.5*π))

initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
storage = simulate!(mech, 3.0, snake_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
# visualize(mech, storage, vis = vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Twister" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

@testset "Energy: Twister" begin
    for jointtype in jointtypes
        # @show jointtype
        mech = get_mechanism(:twister, timestep=timestep0, gravity=gravity0, Nb = Nb0, spring=spring0, damper=damper0,
            jointtype = jointtype, contact = false, r = 0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3] * timestep0
        q10 = UnitQuaternion(RotX(0.5*π))

        initialize!(mech, :twister, q1 = q10, v = v0, ω = ω0)
        storage = simulate!(mech, 3.0, snake_controller!, record = true, verbose = false, opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        ke0 = kinetic_energy(mech, storage)[start0:end]
        pe0 = potential_energy(mech, storage)[start0:end]
        me0 = mechanical_energy(mech, storage)[start0:end]

        # plot([(i-1)*timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
        # plot([(i-1)*timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
        # plot([(i-1)*timestep0 for i in 1:length(me0)], me0 .- me0[1])

        # Test mechanical energy conservation
        @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
        norm((me0 .- me0[1]) ./ mean(me0), Inf)
    end
end
