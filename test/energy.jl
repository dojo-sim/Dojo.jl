# Parameters
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

# Controllers
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

################################################################################
# DICE
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
@testset "Dice" begin
    gravity0 = -10.0

    mech = get_mechanism(:block,
        timestep=timestep0,
        gravity=gravity0,
        contact=false)

    v0 = [1,2,3.0]
    ω0 = [1,1,1.0]
    initialize!(mech, :block,
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-8
end


################################################################################
# SINGLE PENDULUM
################################################################################
# single body
# initial angular velocity
# gravity
# no damper
# no control
################################################################################
@testset "Pendulum" begin
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-2
end

################################################################################
# SLIDER 1
################################################################################
# single body
# no initial linear velocity
# no gravity
# no damper
# no control
################################################################################
@testset "Slider 1" begin
    gravity0 = 0.0
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


    # Test maximum amplitude and velocity
    @test norm(maximum([x[3] for x in storage.x[1]]) - zmax + 0.5) < 1.0e-4
    @test norm(maximum([vl[3] for vl in storage.vl[1]]) - vmax) < 1.0e-4

    # Test mechanical energy conservation
    # With spring the error is or the order of timestep, but there is no drift
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-3
    # @show norm((me0 .- me0[1]) ./ mean(me0), Inf)
end

################################################################################
# SLIDER 2
################################################################################
# single body
# no initial linear velocity
# gravity
# no damper
# no control
################################################################################
@testset "Slider 2" begin
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-6
end

################################################################################
# SLIDER 3
################################################################################
# single body
# no initial linear velocity
# gravity
# no damper
# no control
################################################################################
@testset "Slider 3" begin
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-3
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
@testset "Humanoid" begin
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 2.0e-3
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
@testset "Atlas" begin
    gravity0 = 0.0
    spring0 = 1.0
    damper0 = 0.0
    mech = get_mechanism(:atlas,
        timestep=timestep0,
        gravity=gravity0,
        spring=spring0,
        parse_damper=false,
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 3.0e-3
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
@testset "Quadruped" begin
    gravity0 = 0.0
    spring0 = 1.0
    damper0 = 0.0

    mech = get_mechanism(:quadruped,
        timestep=timestep0,
        gravity=gravity0,
        spring=spring0,
        parse_damper=false,
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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-2
end

################################################################################
# 5-LINK SNAKE
################################################################################
# multiple bodies
# initial linear and angular velocity
# no gravity
# with spring and damper
# with control
################################################################################
@testset "Snake" begin
    Nb0 = 5
    gravity0 = 0.0
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
    q10 = Dojo.RotX(0.5*π)

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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-3


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
        q10 = Dojo.RotX(0.5*π)

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
        @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-2
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
@testset "Twister" begin
    Nb0 = 5
    gravity0 = 0.0
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
    q10 = Dojo.RotX(0.5*π)

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
    @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-3


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
        q10 = Dojo.RotX(0.5*π)

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
        @test norm((me0 .- me0[1]) ./ mean(me0), Inf) < 1.0e-2
        norm((me0 .- me0[1]) ./ mean(me0), Inf)
    end
end
