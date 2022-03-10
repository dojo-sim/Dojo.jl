@testset "FiniteDiff comparison" begin
    q = rand(QuatRotation).q
    @test norm(Dojo.dmrpdq(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.mrp, Dojo.vector(q)), Inf) < 1.0e-5
    @test norm(Dojo.daxisdq(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.axis, Dojo.vector(q))) < 1.0e-5
    @test norm(Dojo.drotation_vectordq(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.rotation_vector, Dojo.vector(q))) < 1.0e-5

    q = one(Quaternion{Float64})
    @test norm(Dojo.dmrpdq(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.mrp, Dojo.vector(q)), Inf) < 1.0e-5
    @test norm(Dojo.drotation_vectordq(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.rotation_vector, Dojo.vector(q))) < 1.0e-5
end
