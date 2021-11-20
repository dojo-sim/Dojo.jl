x0 = srand(3)
v0 = srand(3)
q0 = UnitQuaternion(rand(4)...)
ϕ0 = srand(3)
Δt0 = 0.01
x1 = getx3(x0, v0, Δt0)
q1 = getq3(q0, ϕ0, Δt0)

@testset "Integrator" begin
    @test norm(FiniteDiff.finite_difference_jacobian(x0 -> getx3(x0, v0, Δt0), x0) - Dojo.∂integrator∂x(), Inf) < 1e-8
    @test norm(FiniteDiff.finite_difference_jacobian(v0 -> getx3(x0, v0, Δt0), v0) - Dojo.∂integrator∂v(Δt0), Inf) < 1e-8
    ∂q0 = FiniteDiff.finite_difference_jacobian(q0 -> vector(getq3(UnitQuaternion(q0..., false), ϕ0, Δt0)), vector(q0))
    ∂q1 = Dojo.∂integrator∂q(q0, ϕ0, Δt0, attjac = false)
    @test norm(∂q0 - ∂q1, Inf) < 1e-8
    @test norm(FiniteDiff.finite_difference_jacobian(ϕ0 -> vector(getq3(q0, ϕ0, Δt0)), ϕ0) - Dojo.∂integrator∂ϕ(q0, ϕ0, Δt0), Inf) < 1e-8
end
