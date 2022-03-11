# Controller
function controller!(mechanism, k; 
    U=0.5, 
    timestep=0.01)
    for joint in mechanism.joints
        nu = input_dimension(joint)
        u = (nu <= 5 && k ∈ (1:100)) * U * timestep * sones(nu)
        set_input!(joint, u)
    end
    return
end
nocontrol!(mechanism, k) = controller!(mechanism, k, U=0.0)

# Parameters
ϵ0 = 1e-12
timestep0 = 0.01
gravity0= 0.0
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

################################################################################
# Box
################################################################################
# single body
# initial linear and angular velocity
# no gravity
# no spring and damper
# no control
################################################################################
@testset "Box" begin
    mech = get_mechanism(:block, 
        timestep=timestep0, 
        gravity=gravity0, 
        contact=false)

    v0 = [1,2,3.0]
    ω0 = [10,10,10.0]
    initialize!(mech, :block, 
        v=v0, ω=ω0)

    storage = simulate!(mech, 5.0, nocontrol!, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[5:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-7)
    @test all(norm.(mang0, Inf) .< 1.0e-7)
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
@testset "Pendulum" begin
    mech = get_mechanism(:pendulum, 
        timestep=timestep0, 
        gravity=gravity0)

    ϕ0 = 0.7
    ω0 = 5.0
    initialize!(mech, :pendulum, 
        angle=ϕ0, 
        angular_velocity=ω0)

    storage = simulate!(mech, 5.0, nocontrol!, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[10:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    # @test all(norm.(mlin0, Inf) .< 1e-11)
    @test all(norm.(mang0, Inf) .< 1.0e-7)
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
    spring0 = 1.0
    damper0 = 1.0
    mech = get_mechanism(:humanoid, 
        timestep=timestep0, 
        gravity=gravity0, 
        spring=spring0, 
        damper=damper0, 
        contact_feet=false,
        contact_body=false)

    initialize!(mech, :humanoid)
    bodies = mech.bodies
    set_maximal_velocities!.(bodies, ω=1e-0rand(3))

    storage = simulate!(mech, 10.0, controller!, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[1:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-8)
    @test all(norm.(mang0, Inf) .< 1.0e-8)
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
    spring0 = 10.0
    damper0 = 1.0
    mech = get_mechanism(:atlas, 
        timestep=timestep0, 
        gravity=gravity0, 
        spring=spring0, 
        damper=damper0, 
        contact_feet=false,
        contact_body=false)

    initialize!(mech, :atlas)
    storage = simulate!(mech, 5.0, controller!, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[1:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-8)
    @test all(norm.(mang0, Inf) .< 1.0e-8)
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
    spring0 = 0.3
    damper0 = 0.1
    mech = get_mechanism(:quadruped, 
        timestep=timestep0, 
        gravity=gravity0, 
        spring=spring0, 
        damper=damper0, 
        contact_feet=false,
        contact_body=false)

    initialize!(mech, :quadruped)
    storage = simulate!(mech, 5.0, controller!, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[1:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-8)
    @test all(norm.(mang0, Inf) .< 1.0e-8)
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
@testset "Snake" begin
    Nb0 = 5
    spring0 = 0.0 * 4e0
    damper0 = 0.0 * 2e+1

    mech = get_mechanism(:snake, 
        timestep=timestep0, 
        gravity=gravity0, 
        num_bodies=Nb0, 
        spring=spring0, 
        damper=damper0,
        joint_type=:Revolute, 
        contact=false, 
        radius=0.05);

    v0 = 100.0 * [1, 2, 3] * timestep0
    ω0 = 100.0 * [1, 2, 3.0] * timestep0
    q10 = Dojo.RotX(0.5*π)

    initialize!(mech, :snake, 
        base_orientation=q10, 
        base_linear_velocity=v0, 
        base_angular_velocity=ω0)
    storage = simulate!(mech, 1.50, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[5:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-8)
    @test all(norm.(mang0, Inf) .< 1.0e-8)

    
    for joint_type in joint_types
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
        ω0 = 10.0 * [1, 2, 3.0] * timestep0
        q10 = Dojo.RotX(0.5 * π)

        initialize!(mech, :snake, 
            base_orientation=q10, 
            base_linear_velocity=v0, 
            base_angular_velocity=ω0)
        storage = simulate!(mech, 1.50, controller!, 
            record=true, 
            verbose=false, 
            opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        # @show mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
        # @show mang0= [Vector(m - m0[1])[4:6] for m in m0]
        # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
        # plt = plot()
        # plot!([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
        # display(plt)
        @test all(norm.(mlin0, Inf) .< 1.0e-8)
        @test all(norm.(mang0, Inf) .< 1.0e-8)
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
@testset "Twister" begin
    Nb0 = 5
    spring0 = 1.0 * 4e0
    damper0 = 1.0 * 2e+1

    mech = get_mechanism(:twister, 
        timestep=timestep0, 
        gravity=gravity0, 
        num_bodies=Nb0, 
        spring=spring0, 
        damper=damper0,
        joint_type=:FixedOrientation, 
        contact=false, 
        radius=0.05);

    v0 = 100.0 * [1, 2, 3] * timestep0
    ω0 = 100.0 * [1, 2, 3.0] * timestep0
    q10 = Dojo.RotX(0.5*π)

    initialize!(mech, :twister, 
        base_orientation=q10, 
        base_linear_velocity=v0, 
        base_angular_velocity=ω0)
    storage = simulate!(mech, 2.5, 
        record=true, 
        verbose=false, 
        opts=SolverOptions(rtol=ϵ0, btol=ϵ0))

    # visualize(mech, storage, vis = vis)

    m0 = momentum(mech, storage)[5:end]
    mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
    mang0= [Vector(m - m0[1])[4:6] for m in m0]
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mlin0...)')
    # plot([(i-1) * timestep0 for i in 1:length(m0)], hcat(mang0...)')
    @test all(norm.(mlin0, Inf) .< 1.0e-8)
    @test all(norm.(mang0, Inf) .< 1.0e-8)

    
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
        ω0 = 10.0 * [1, 2, 3.0] * timestep0
        q10 = Dojo.RotX(0.5*π)
        initialize!(mech, :twister, 
            base_orientation=q10, 
            base_linear_velocity=v0, 
            base_angular_velocity=ω0)
        storage = simulate!(mech, 1.50, 
            record=true, 
            verbose=false, 
            opts=SolverOptions(rtol=ϵ0, btol=ϵ0))
        # visualize(mech, storage, vis = vis)

        m0 = momentum(mech, storage)[5:end]
        mlin0 = [Vector(m - m0[1])[1:3] for m in m0]
        mang0= [Vector(m - m0[1])[4:6] for m in m0]

        @test all(norm.(mlin0, Inf) .< 1.0e-8)
        @test all(norm.(mang0, Inf) .< 1.0e-8)
    end
end