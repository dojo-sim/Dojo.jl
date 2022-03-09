# Data
ϵ0 = 1.0e-12
timestep0 = 1.0e-2
start0 = Int(floor(1 / timestep0)) + 1

joint_types = [
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
function controller!(mechanism, k; 
    U=0.5, 
    timestep=timestep0)

    N = Int(floor(1 / timestep))
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U=0.0)

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

mech = get_mechanism(:box, 
    timestep=timestep0, 
    gravity=gravity0, 
    contact=false)

v0 = [1,2,3.0]
ω0 = [1,1,1.0]
initialize!(mech, :box, 
    v=v0, 
    ω=ω0)

storage = simulate!(mech, 5.0, nocontrol!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
# With spring the error is or the order of timestep, but there is no drift
@testset "Energy: Dice" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-8
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
gravity0 = 0.0
spring0 = 1.0
damper0 = 0.0

mech = get_mechanism(:pendulum, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0, 
    damper=damper0)

ϕ0 = 0.5π
ω0 = 0.0

initialize!(mech, :pendulum, 
    angle=ϕ0, 
    angular_velocity=ω0)
storage = simulate!(mech, 25.0, controller!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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
mech = get_mechanism(:slider, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0, 
    damper=damper0)

z0 = 0.5
initialize!(mech, :slider, 
    position=z0)

# Analytical
pbody = mech.bodies[1]
zmax = z0
vmax = z0 * sqrt(spring0 / pbody.mass)
pe_max = 0.5 * spring0 * zmax^2
ke_max = 0.5 * pbody.mass * vmax^2

storage = simulate!(mech, 5.0,  nocontrol!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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

mech = get_mechanism(:slider, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0, 
    damper=damper0)

z0 = 0.5
initialize!(mech, :slider, 
    position=z0)

storage = simulate!(mech, 1.5,  nocontrol!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
# For gravity the conservation is perfect
@testset "Energy: Slider" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-6
end
norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-6


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
spring0 = 1.0
damper0 = 0.0

mech = get_mechanism(:slider, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0, 
    damper=damper0)

z0 = 0.1
initialize!(mech, :slider, 
    position=z0)

storage = simulate!(mech, 10.0,  nocontrol!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[start0])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[start0])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[start0])

# Test mechanical energy conservation
@testset "Energy: Slider" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
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
function humanoid_controller!(mechanism, k; 
    U=0.05, 
    timestep=timestep0)

    N = Int(floor(1 / timestep))
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

gravity0 = 0.0
spring0 = 1.0
damper0 = 0.0

mech = get_mechanism(:humanoid, 
    timestep=timestep0, 
    gravity=gravity0,
    spring=spring0, 
    damper=damper0, 
    contact_feet=false,
    contact_body=false)

initialize!(mech, :humanoid)
bodies = mech.bodies
for body in mech.bodies
    set_maximal_velocities!(body, 
        ω=0.5 * rand(3))
    # set_maximal_velocities!(body, v=1.0 * rand(3))
end

for joint in mech.joints 
    joint.rotational.spring_type = :linear 
end

storage = simulate!(mech, 3.0, humanoid_controller!, 
    record=true,
    verbose=false, 
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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
mech = get_mechanism(:atlas, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0,
    damper=damper0, 
    contact_feet=false,
    contact_body=false)

initialize!(mech, :atlas)
bodies = mech.bodies
set_maximal_velocities!.(bodies, 
    ω=1.0 * rand(3))

storage = simulate!(mech, 5.0, humanoid_controller!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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
function quadruped_controller!(mechanism, k; 
    U=0.05, 
    timestep=timestep0)

    N = Int(floor(1.0 / timestep))
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

spring0 = 1.0
damper0 = 0.0

mech = get_mechanism(:quadruped, 
    timestep=timestep0, 
    gravity=gravity0, 
    spring=spring0,
    damper=damper0, 
    contact_feet=false, 
    contact_body=false,
    limits=false)

initialize!(mech, :quadruped)
storage = simulate!(mech, 5.0, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
    
# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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
function snake_controller!(mechanism, k; 
    U=0.05, 
    timestep=timestep0)
    N = Int(floor(1.0 / timestep))
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

Nb0 = 5
spring0 = 0.01
damper0 = 0.0

mech = get_mechanism(:snake, 
    timestep=timestep0, 
    gravity=gravity0, 
    num_bodies=Nb0,
    spring=spring0, 
    damper=damper0, 
    joint_type=:Revolute, 
    contact=false, 
    radius=0.05);

v0 = 10.0 * [1, 2, 3] * timestep0
ω0 = 10.0 * [1, 2, 3] * timestep0
q10 = Quaternion(RotX(0.5*π))

initialize!(mech, :snake, 
    base_orientation=q10, 
    base_linear_velocity=v0, 
    base_angular_velocity=ω0)
storage = simulate!(mech, 3.0, snake_controller!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Snake" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

@testset "Energy: Snake" begin
    for joint_type in joint_types
        # @show joint_type
        mech = get_mechanism(:snake, 
            timestep=timestep0, 
            gravity=gravity0, 
            num_bodies=Nb0,
            spring=spring0, 
            damper=damper0, 
            joint_type=joint_type, 
            contact=false, 
            radius=0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3] * timestep0
        q10 = Quaternion(RotX(0.5*π))

        initialize!(mech, :snake, 
            base_orientation=q10, 
            base_linear_velocity=v0, 
            base_angular_velocity=ω0)
        storage = simulate!(mech, 3.0, snake_controller!, 
            record=true, 
            verbose=false,
            opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

        # visualize(mech, storage, vis=vis)

        ke0 = kinetic_energy(mech, storage)[start0:end]
        pe0 = potential_energy(mech, storage)[start0:end]
        me0 = mechanical_energy(mech, storage)[start0:end]

        # plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
        # plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
        # plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

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
function twister_controller!(mechanism, k; 
    U=0.01, 
    timestep=timestep0)
    N = Int(floor(1.0 / timestep))
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:N)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end

Nb0 = 5
spring0 = 0.01
damper0 = 0.0

mech = get_mechanism(:twister, 
    timestep=timestep0, 
    gravity=gravity0, 
    num_bodies=Nb0,
    spring=spring0, 
    damper=damper0, 
    joint_type=:Revolute, 
    contact=false, 
    radius=0.05);

v0 = 10.0 * [1, 2, 3] * timestep0
ω0 = 10.0 * [1, 2, 3] * timestep0
q10 = Quaternion(RotX(0.5*π))

initialize!(mech, :twister, 
    base_orientation=q10, 
    base_linear_velocity=v0, 
    base_angular_velocity=ω0)
storage = simulate!(mech, 3.0, twister_controller!, 
    record=true, 
    verbose=false,
    opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

# visualize(mech, storage, vis=vis)

ke0 = kinetic_energy(mech, storage)[start0:end]
pe0 = potential_energy(mech, storage)[start0:end]
me0 = mechanical_energy(mech, storage)[start0:end]

# plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
# plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
# plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

# Test mechanical energy conservation
@testset "Energy: Twister" begin
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-3
end
norm((me0 .- me0[1]) ./ mean(me0), Inf)

@testset "Energy: Twister" begin
    for joint_type in joint_types
        mech = get_mechanism(:twister, 
            timestep=timestep0, 
            gravity=gravity0, 
            num_bodies=Nb0,
            spring=spring0, 
            damper=damper0, 
            joint_type=joint_type, 
            contact=false, 
            radius=0.05)

        v0 = 10.0 * [1, 2, 3] * timestep0
        ω0 = 10.0 * [1, 2, 3] * timestep0
        q10 = Quaternion(RotX(0.5*π))

        initialize!(mech, :twister, 
            base_orientation=q10, 
            base_linear_velocity=v0, 
            base_angular_velocity=ω0)
        storage = simulate!(mech, 3.0, snake_controller!, 
            record=true, 
            verbose=false,
            opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

        # visualize(mech, storage, vis=vis)

        ke0 = kinetic_energy(mech, storage)[start0:end]
        pe0 = potential_energy(mech, storage)[start0:end]
        me0 = mechanical_energy(mech, storage)[start0:end]

        # plot([(i-1) * timestep0 for i in 1:length(ke0)], ke0 .- ke0[1])
        # plot([(i-1) * timestep0 for i in 1:length(pe0)], pe0 .- pe0[1])
        # plot([(i-1) * timestep0 for i in 1:length(me0)], me0 .- me0[1])

        # Test mechanical energy conservation
        @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1e-2
        norm((me0 .- me0[1]) ./ mean(me0), Inf)
    end
end
