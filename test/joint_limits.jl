@testset "Pendulum" begin
    mech = get_mechanism(:pendulum;
        timestep=0.01, 
        gravity=-9.81, 
        springs=0.0, 
        dampers=0.0,
        joint_limits=Dict([
        (:joint, [0.25 * π, 1.0 * π]),])
    )

    initialize!(mech, :pendulum; 
        angle=0.4 * π)
    storage = simulate!(mech, 1.0; 
        record=true, 
        opts=SolverOptions(btol=1.0e-5, rtol=1.0e-5))

    @test norm(Dojo.get_minimal_state(mech)[1] - 0.25 * π) < 1.0e-3
end
