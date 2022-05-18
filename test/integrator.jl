@testset "ForwardDiff comparison" begin
    Random.seed!(100)
    x0 = srand(3)
    v0 = srand(3)
    q0 = rand(QuatRotation).q
    ω0 = srand(3)
    timestep0 = 0.01
    x1 = Dojo.next_position(x0, v0, timestep0)
    q1 = Dojo.next_orientation(q0, ω0, timestep0)

    # ∇v
    FD∂i∂v = ForwardDiff.jacobian(v0 -> Dojo.next_position(x0, v0, timestep0), v0)
    @test norm(FD∂i∂v - Dojo.linear_integrator_jacobian_velocity(x0, v0, timestep0), Inf) < 1.0e-8
    # ∇ω
    FD∂i∂ω = ForwardDiff.jacobian(ω0 -> Dojo.vector(Dojo.next_orientation(q0, ω0, timestep0)), ω0)
    @test norm(FD∂i∂ω - Dojo.rotational_integrator_jacobian_velocity(q0, ω0, timestep0), Inf) < 1.0e-8

    # ∇x
    FD∂i∂x = ForwardDiff.jacobian(x0 -> Dojo.next_position(x0, v0, timestep0), x0)
    @test norm(FD∂i∂x - Dojo.linear_integrator_jacobian_position(x0, v0, timestep0), Inf) < 1.0e-8
    # ∇q
    FD∂i∂q = ForwardDiff.jacobian(q0 ->
        Dojo.vector(Dojo.next_orientation(Quaternion(q0...), ω0, timestep0)),
        Dojo.vector(q0))
    @test norm(FD∂i∂q - Dojo.rotational_integrator_jacobian_orientation(q0, ω0, timestep0, attjac=false), Inf) < 1.0e-8

    # ∇vω
    vel0 = [v0; ω0]
    FDvel = ForwardDiff.jacobian(vel0 ->
        [Dojo.next_position(x0, SVector{3}(vel0[1:3]), timestep0);
        Dojo.vector(Dojo.next_orientation(q0, SVector{3}(vel0[4:6]), timestep0))], vel0)
    @test norm(FDvel - Dojo.integrator_jacobian_velocity(x0, v0, q0, ω0, timestep0), Inf) < 1.0e-8

    # ∇xq
    con0 = [x0; Dojo.vector(q0)]
    FDcon = ForwardDiff.jacobian(con0 ->
        [Dojo.next_position(SVector{3}(con0[1:3]), v0, timestep0);
        Dojo.vector(Dojo.next_orientation(Quaternion(con0[4:7]...), ω0, timestep0))], con0)
    @test norm(FDcon - Dojo.integrator_jacobian_configuration(x0, v0, q0, ω0, timestep0, attjac=false), Inf) < 1.0e-8
    @test norm(FDcon * cat(I(3), Dojo.LVᵀmat(q0), dims=(1,2)) - Dojo.integrator_jacobian_configuration(x0, v0, q0, ω0, timestep0, attjac=true), Inf) < 1.0e-8
end
