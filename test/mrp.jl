@testset "FiniteDiff comparison" begin
    q = [1,2,3,4.0]
    q = Quaternion(q ./ norm(q)...,true)
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


@testset "axes pair to quaternion" begin
    n1 = [0,0,0.0]
    n2 = [0,0,0.0]
    q = Dojo.axes_pair_to_quaternion(n1, n2)
    @test norm(n2 - Dojo.vector_rotate(n1, q), Inf) < 1e-10

    n1 = rand(3)
    n1 /= norm(n1)
    n2 = n1
    q = Dojo.axes_pair_to_quaternion(n1, n2)
    @test norm(n2 - Dojo.vector_rotate(n1, q), Inf) < 1e-10

    n1 = rand(3)
    n1 /= norm(n1)
    n2 = -n1
    q = Dojo.axes_pair_to_quaternion(n1, n2)
    @test norm(n2 - Dojo.vector_rotate(n1, q), Inf) < 1e-5

    n1 = rand(3)
    n1 /= norm(n1)
    n2 = rand(3)
    n2 /= norm(n2)
    q = Dojo.axes_pair_to_quaternion(n1, n2)
    @test norm(n2 - Dojo.vector_rotate(n1, q), Inf) < 1e-10
end
