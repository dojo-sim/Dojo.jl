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

@testset "rotation matrix jacobian" begin
    p0 = rand(3)
    q0 = rand(QuatRotation).q
    J0 = Dojo.∂rotation_matrix∂q(q0, p0)
    J2 = Dojo.∂rotation_matrix∂q(q0, p0, attjac=true)
    J1 = FiniteDiff.finite_difference_jacobian(q0 -> Dojo.rotation_matrix(Quaternion(q0...,true))*p0, Dojo.vector(q0))
    J3 = J1 * Dojo.LVᵀmat(q0)
    @test norm(J0 - J1, Inf) < 1e-6
    @test norm(J2 - J3, Inf) < 1e-6

    J4 = Dojo.∂rotation_matrix_inv∂q(q0, p0)
    J6 = Dojo.∂rotation_matrix_inv∂q(q0, p0, attjac=true)
    J5 = FiniteDiff.finite_difference_jacobian(q0 -> Dojo.rotation_matrix(inv(Quaternion(q0...,true)))*p0, Dojo.vector(q0))
    J7 = J5 * Dojo.LVᵀmat(q0)
    @test norm(J4 - J5, Inf) < 1e-6
    @test norm(J6 - J7, Inf) < 1e-6
end
