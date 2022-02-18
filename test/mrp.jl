@testset "Modified Rodrigues Parameters" begin
    q = rand(UnitQuaternion)
    @test norm(Dojo.∂mrp∂q(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.mrp, Dojo.vector(q)), Inf) < 1.0e-6
    @test norm(Dojo.∂axis∂q(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.axis, Dojo.vector(q))) < 1.0e-6
    @test norm(Dojo.∂rotation_vector∂q(Dojo.vector(q)) -
        FiniteDiff.finite_difference_jacobian(Dojo.rotation_vector, Dojo.vector(q))) < 1.0e-6
    #TODO: add tests for q = [1.0, 0.0, 0.0, 0.0]
end
