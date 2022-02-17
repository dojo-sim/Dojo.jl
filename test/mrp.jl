@testset "Modified Rodrigues Parameters" begin
    q = rand(UnitQuaternion)
    @test norm(∂mrp∂q(vector(q)) - FiniteDiff.finite_difference_jacobian(mrp, vector(q)), Inf) < 1.0e-6
    @test norm(∂axis∂q(vector(q)) - FiniteDiff.finite_difference_jacobian(axis, vector(q))) < 1.0e-6
    @test norm(∂rotation_vector∂q(vector(q)) - FiniteDiff.finite_difference_jacobian(rotation_vector, vector(q))) < 1.0e-6
    #TODO: add tests for q = [1.0, 0.0, 0.0, 0.0]
end
