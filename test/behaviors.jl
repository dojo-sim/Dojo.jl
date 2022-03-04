@testset "Behavior: Quadruped simulation" begin
    mech = get_mechanism(:quadruped,
        timestep=0.05,
        gravity=-9.81,
        friction_coefficient=0.8,
        damper=1000.0,
        spring=30.0)

    initialize!(mech, :quadruped)

    try
        storage = simulate!(mech, 5.0,
            record=true,
            verbose=false)
        @test true
    catch
        @test false
    end
end

@testset "Behavior: Box toss" begin
    for timestep in [0.10, 0.05, 0.01, 0.005]
        mech = get_mechanism(:box,
            timestep=timestep,
            gravity=-9.81,
            friction_coefficient = 0.1)

        initialize!(mech, :box,
            x=[0.0, 0.0, 0.5],
            v=[1.0, 1.5, 1.0],
            ω=[5.0, 4.0, 2.0] .* timestep)
        storage = simulate!(mech, 5.0,
            record=true,
            opts=SolverOptions(btol=1e-6, rtol=1e-6,
            verbose=false))

        @test norm(storage.v[1][end], Inf) < 1e-12
        @test norm(storage.x[1][end][3] - 0.25, Inf) < 1e-3
    end
end

@testset "Behavior: Four-bar linkage" begin
    for timestep in [0.10, 0.05, 0.01, 0.005]
        mech = Dojo.get_mechanism(:fourbar,
            model="fourbar",
            timestep=timestep)
        Dojo.initialize!(mech, :fourbar,
            angle=0.25)
        loopjoints = mech.joints[end:end]
        Dojo.root_to_leaves_ordering(mech) == [2, 7, 3, 6, 1, 8, 4, 9]

        # Simulation
        function ctrl!(m, t)
            Dojo.set_input!(m, 1.0 * m.timestep * SVector(rand(), -rand(), 0.0, 0.0, 0.0))
            return nothing
        end
        storage = Dojo.simulate!(mech, 5.0, ctrl!, verbose=false, record=true)

        min_coords = Dojo.get_minimal_coordinates(mech)
        @test norm(min_coords[5] - +min_coords[4], Inf) < 1.0e-5
        @test norm(min_coords[5] - -min_coords[3], Inf) < 1.0e-5
        @test norm(min_coords[5] - (min_coords[2] - min_coords[1]), Inf) < 1.0e-5
    end
end

@testset "Behavior: Tennis Racket" begin
    # Simulation
    timestep=0.01
    gravity=0.0
    mech = Dojo.get_mechanism(:dzhanibekov,
            timestep=timestep,
            gravity=gravity);

    # Simulate
    Dojo.initialize!(mech, :dzhanibekov,
        angular_velocity=[15.0; 0.01; 0.0])
    storage = Dojo.simulate!(mech, 4.0,
        record=true,
        verbose=false)

    # The x position of the side body oscillate between a positive and negative
    # value when we observe the Dhanibekov effect. Otherwise, it always remains in the positive position.
    @test Dojo.minimum([x[1] for x in storage.x[2]]) < -0.05
end
