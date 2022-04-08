using Test
using FiniteDiff
using LinearAlgebra

@testset "coulomb direction" begin
    v = srand(3)
    Dojo.coulomb_direction(v)
    J0 = Dojo.∂coulomb_direction∂v(v)
    J1 = FiniteDiff.finite_difference_jacobian(v-> Dojo.coulomb_direction(v), v)
    @test norm(J0 - J1, Inf) < 1e-5
end

# solmat


# datamat
