@testset "Integrator" begin
    Random.seed!(100)
    x0 = srand(3)
    v0 = srand(3)
    q0 = UnitQuaternion(rand(4)...)
    ϕ0 = srand(3)
    timestep0 = 0.01
    x1 = Dojo.next_position(x0, v0, timestep0)
    q1 = Dojo.next_orientation(q0, ϕ0, timestep0)

    @test norm(FiniteDiff.finite_difference_jacobian(x0 -> Dojo.next_position(x0, v0, timestep0), x0) - Dojo.∂integrator∂x(), Inf) < 1e-8
    @test norm(FiniteDiff.finite_difference_jacobian(v0 -> Dojo.next_position(x0, v0, timestep0), v0) - Dojo.∂integrator∂v(timestep0), Inf) < 1e-8
    ∂q0 = FiniteDiff.finite_difference_jacobian(q0 -> Dojo.vector(Dojo.next_orientation(UnitQuaternion(q0..., false), ϕ0, timestep0)), Dojo.vector(q0))
    ∂q1 = Dojo.∂integrator∂q(q0, ϕ0, timestep0, attjac = false)
    @test norm(∂q0 - ∂q1, Inf) < 1e-8
    @test norm(FiniteDiff.finite_difference_jacobian(ϕ0 -> Dojo.vector(Dojo.next_orientation(q0, ϕ0, timestep0)), ϕ0) - Dojo.∂integrator∂ϕ(q0, ϕ0, timestep0), Inf) < 1e-8
end
