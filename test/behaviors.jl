@testset "Quadruped" begin
    mech = get_mechanism(:quadruped;
        timestep=0.05,
        gravity=-9.81,
        friction_coefficient=0.8,
        parse_springs=false,
        parse_dampers=false,
        dampers=1000.0,
        springs=30.0)

    initialize!(mech, :quadruped)

    storage = simulate!(mech, 5.0,
            record=true,
            verbose=false)

    res = get_sdf(mech, storage) # distance from floor to each contact
    @test minimum(minimum([min.(0.0, r) for r in res])) >= 0.0
end

@testset "Box toss" begin
    for timestep in [0.10, 0.05, 0.01, 0.005]
        mech = get_mechanism(:block;
            timestep,
            gravity=-9.81,
            friction_coefficient = 0.1)

        initialize!(mech, :block,
            position=[0.0, 0.0, 0.5],
            velocity=[1.0, 1.5, 1.0],
            angular_velocity=[5.0, 4.0, 2.0] .* timestep)
        storage = simulate!(mech, 5.0,
            record=true,
            opts=SolverOptions(btol=1e-6, rtol=1e-6,
            verbose=false))

        @test norm(storage.v[1][end], Inf) < 1.0e-8
        @test norm(storage.x[1][end][3] - 0.25, Inf) < 1.0e-3
    end
end

@testset "Box external force" begin
    mech = get_mechanism(:block;gravity=0, contact=false, mass=1)
    mech.bodies[1].inertia = diagm(ones(3))

    initialize!(mech, :block; position=zeros(3), orientation=one(Quaternion), velocity=zeros(3), angular_velocity=zeros(3))
    set_external_force!(mech.bodies[1];force=[1;0;0])
    storage = simulate!(mech, 1.0, record=true)
    @test norm(mech.bodies[1].state.vsol[1][1] - 1.0) < 1.0e-3

    initialize!(mech, :block; position=zeros(3), orientation=one(Quaternion), velocity=zeros(3), angular_velocity=zeros(3))
    set_external_force!(mech.bodies[1];torque=[1;0;0])
    storage = simulate!(mech, 1.0, record=true)
    @test norm(mech.bodies[1].state.ωsol[1][1] - 1.0) < 1.0e-3
end

@testset "Four-bar linkage" begin
    for timestep in [0.10, 0.05, 0.01, 0.005]
        mech = get_mechanism(:fourbar;
            timestep,
            parse_springs=false,
			parse_dampers=false,
            dampers=0.0)
        Dojo.initialize!(mech, :fourbar,
            inner_angle=0.25)
        loopjoints = mech.joints[end:end]
        Dojo.root_to_leaves_ordering(mech) == [2, 7, 3, 6, 1, 8, 4, 9]

        # Simulation
        function ctrl!(m, t)
            Dojo.set_input!(m, 1.0 * SVector(rand(), -rand(), 0.0, 0.0, 0.0))
            return nothing
        end
        storage = Dojo.simulate!(mech, 5.0, ctrl!, verbose=false, record=true)

        min_coords = Dojo.get_minimal_coordinates(mech)
        @test norm(min_coords[5] - +min_coords[4], Inf) < 1.0e-5
        @test norm(min_coords[5] - -min_coords[3], Inf) < 1.0e-5
        @test norm(min_coords[5] - (min_coords[2] - min_coords[1]), Inf) < 1.0e-5
    end
end

@testset "Dzhanibekov" begin
    # Simulation
    timestep=0.01
    gravity=0.0
    mech = get_mechanism(:dzhanibekov;
            timestep,
            gravity);

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
