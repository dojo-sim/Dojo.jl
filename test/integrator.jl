@testset "Integrator" begin
    Random.seed!(100)
    x0 = srand(3)
    v0 = srand(3)
    q0 = rand(QuatRotation).q
    ϕ0 = srand(3)
    timestep0 = 0.01
    x1 = Dojo.next_position(x0, v0, timestep0)
    q1 = Dojo.next_orientation(q0, ϕ0, timestep0)

    # ∇v
    FD∂i∂v = FiniteDiff.finite_difference_jacobian(v0 -> Dojo.next_position(x0, v0, timestep0), v0)
    @test norm(FD∂i∂v - Dojo.linear_integrator_jacobian_velocity(timestep0), Inf) < 1e-8
    # ∇ϕ
    FD∂i∂ϕ = FiniteDiff.finite_difference_jacobian(ϕ0 -> Dojo.vector(Dojo.next_orientation(q0, ϕ0, timestep0)), ϕ0)
    @test norm(FD∂i∂ϕ - Dojo.rotational_integrator_jacobian_velocity(q0, ϕ0, timestep0), Inf) < 1e-8

    # ∇x
    FD∂i∂x = FiniteDiff.finite_difference_jacobian(x0 -> Dojo.next_position(x0, v0, timestep0), x0)
    @test norm(FD∂i∂x - Dojo.linear_integrator_jacobian_position(), Inf) < 1e-8
    # ∇q
    FD∂i∂q = FiniteDiff.finite_difference_jacobian(q0 ->
        Dojo.vector(Dojo.next_orientation(Quaternion(q0...), ϕ0, timestep0)),
        Dojo.vector(q0))
    @test norm(FD∂i∂q - Dojo.rotational_integrator_jacobian_orientation(q0, ϕ0, timestep0, attjac=false), Inf) < 1e-8

    # ∇vϕ
    vel0 = [v0; ϕ0]
    FDvel = FiniteDiff.finite_difference_jacobian(vel0 ->
        [Dojo.next_position(x0, SVector{3}(vel0[1:3]), timestep0);
        Dojo.vector(Dojo.next_orientation(q0, SVector{3}(vel0[4:6]), timestep0))], vel0)
    @test norm(FDvel - Dojo.integrator_jacobian_velocity(q0, ϕ0, timestep0), Inf) < 1e-8

    # ∇xq
    con0 = [x0; Dojo.vector(q0)]
    FDcon = FiniteDiff.finite_difference_jacobian(con0 ->
        [Dojo.next_position(SVector{3}(con0[1:3]), v0, timestep0);
        Dojo.vector(Dojo.next_orientation(Quaternion(con0[4:7]...), ϕ0, timestep0))], con0)
    @test norm(FDcon - Dojo.integrator_jacobian_configuration(q0, ϕ0, timestep0, attjac=false), Inf) < 1e-8
    @test norm(FDcon * cat(I(3), Dojo.LVᵀmat(q0), dims=(1,2)) - Dojo.integrator_jacobian_configuration(q0, ϕ0, timestep0, attjac=true), Inf) < 1e-8
end
