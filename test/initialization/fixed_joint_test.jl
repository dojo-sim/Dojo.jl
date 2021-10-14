using ConstrainedDynamics
using ConstrainedDynamics: vrotate
using Rotations
using StaticArrays
using LinearAlgebra


Δt = 0.01
length1 = 1.0
width, depth = 1.0, 1.0

for i=1:10
    # Links
    origin = Origin{Float64}()
    link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
    link1.m = 1.0
    link1.J = diagm(ones(3))

    p1 = rand(3)
    p2 = rand(3)
    qoff = rand(UnitQuaternion)


    # Constraints
    joint1 = EqualityConstraint(Fixed(origin, link1; p1=p1, p2=p2, qoffset = qoff))

    links = [link1]
    constraints = [joint1]

    mech = Mechanism(origin, links, constraints, g = 0.)

    setPosition!(mech,joint1,SA_F64[]; iter=false)
    setVelocity!(mech,joint1,SA_F64[])
    # setForce!(mech,joint1,SA_F64[]) # Does not do anything

    storage = simulate!(mech, 10., record = true)

    truex = [(p1 - vrotate(p2,qoff)) for i=0:999]
    trueq = [qoff for i=0:999]

    @test isapprox(norm(minimalCoordinates(mech, joint1) - SA_F64[]), 0.0; atol = 1e-8)
    @test isapprox(sum(norm.(storage.x[1].-truex))/1000, 0.0; atol = 1e-8)
    @test isapprox(sum(norm.(storage.q[1].-trueq))/1000, 0.0; atol = 1e-8)
end