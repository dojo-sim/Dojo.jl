# @testset "FiniteDiff comparison" begin
    q = rand(QuatRotation).q
    @test norm(Dojo.dmrpdq(Dojo.vector(q)) -
        ForwardDiff.jacobian(Dojo.mrp, Dojo.vector(q)), Inf) < 1.0e-5
    @test norm(Dojo.daxisdq(Dojo.vector(q)) -
        ForwardDiff.jacobian(Dojo.axis, Dojo.vector(q))) < 1.0e-5
    @test norm(Dojo.drotation_vectordq(Dojo.vector(q)) -
        ForwardDiff.jacobian(Dojo.rotation_vector, Dojo.vector(q))) < 1.0e-5
    q = one(Quaternion{Float64})
    @test norm(Dojo.dmrpdq(Dojo.vector(q)) -
        ForwardDiff.jacobian(Dojo.mrp, Dojo.vector(q)), Inf) < 1.0e-5
    @test norm(Dojo.drotation_vectordq(Dojo.vector(q)) -
        SA[
            0.0  2.0  0.0  0.0;
            0.0  0.0  2.0  0.0;
            0.0  0.0  0.0  2.0;
        ]) < 1.0e-5
    # TODO For zero rotation, ForwardDiff runs into a numeric singularity
    # @test norm(Dojo.drotation_vectordq(Dojo.vector(q)) -
    #     ForwardDiff.jacobian(Dojo.rotation_vector, Dojo.vector(q))) < 1.0e-5
# end
