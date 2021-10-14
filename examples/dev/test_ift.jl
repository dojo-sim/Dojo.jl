using Test
using LinearAlgebra

@testset "integration scheme gradient" begin

    # integration scheme gradient
    # test
    Δt = 0.01
    x20 = rand(3)
    q20 = rand(4)
    q20 ./= norm(q20)
    v20 = rand(3)
    ω20 = rand(3)
    x30 = st(x20, v20, Δt)
    q30 = sr(q20, ω20, Δt)
    x31 = getx3(x20, v20, Δt)
    q31 = getq3(q20, ω20, Δt)
    q31 = [q31.w, q31.x, q31.y, q31.z]
    @test norm(x30 - x31) < 1e-10
    @test norm(q30 - q31) < 1e-10

    # test Δt = 0
    Δt = 1e-20
    x20 = rand(3)
    q20 = rand(4)
    q20 ./= norm(q20)
    v20 = rand(3)
    ω20 = rand(3)
    x30 = st(x20, v20, Δt)
    q30 = sr(q20, ω20, Δt)
    @test norm(x20 - x30) < 1e-10
    @test norm(q20 - q30) < 1e-10

    # test v2 = 0; ω2 = 0
    Δt = 0.01
    x20 = rand(3)
    q20 = rand(4)
    q20 ./= norm(q20)
    v20 = zeros(3)
    ω20 = zeros(3)
    x30 = st(x20, v20, Δt)
    q30 = sr(q20, ω20, Δt)
    @test norm(x20 - x30) < 1e-10
    @test norm(q20 - q30) < 1e-10

    # test back and forth
    Δt = 0.01
    x10 = rand(3)
    q10 = rand(4)
    q10 ./= norm(q10)
    v0 = zeros(3)
    ω0 = zeros(3)
    x20 = st(x10, v0, Δt)
    q20 = sr(q10, ω0, Δt)
    x30 = st(x20, -v0, Δt)
    q30 = sr(q20, -ω0, Δt)
    @test norm(x10 - x30) < 1e-10
    @test norm(q10 - q30) < 1e-10

    # test derivatives
    Δt = 0.1
    x20 = rand(3)
    q20 = rand(4)
    q20 ./= norm(q20)
    v20 = rand(3)
    ω20 = rand(3)
    srω0 = ∂sr∂ω(q20, ω20, Δt)
    srω1 = ForwardDiff.jacobian(ω -> sr(q20, ω, Δt), ω20)
    @test norm(srω0 - srω1) < 1e-10

    srq0 = ∂sr∂q(q20, ω20, Δt)
    srq1 = ForwardDiff.jacobian(q -> sr(q, ω20, Δt, normalize = false), q20)
    @test norm(srq0 - srq1) < 1e-10

    stx0 = ∂st∂x(x20, v20, Δt)
    stx1 = ForwardDiff.jacobian(x -> st(x, v20, Δt), x20)
    @test norm(stx0 - stx1) < 1e-10

    stv0 = ∂st∂v(x20, v20, Δt)
    stv1 = ForwardDiff.jacobian(v -> st(x20, v, Δt), v20)
    @test norm(stv0 - stv1) < 1e-10
end


@testset "implicit function theorem" begin


end
