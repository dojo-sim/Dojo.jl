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
    axis = rand(3)
    axis = axis/norm(axis)
    qoff = rand(UnitQuaternion)


    # Constraints
    joint1 = EqualityConstraint(Prismatic(origin, link1, axis; p1=p1, p2=p2, qoffset = qoff))

    links = [link1]
    constraints = [joint1]

    mech = Mechanism(origin, links, constraints, g = 0.)

    xθ = rand(1)
    vω = rand(1)
    Fτ = rand(1)

    setPosition!(mech,joint1,xθ;iter=false)
    setVelocity!(mech,joint1,vω)
    function control!(mechanism, k)
        if k==1
            setForce!(mechanism,joint1,Fτ) 
        end
        return
    end

    storage = simulate!(mech, 10., control!, record = true)

    @test isapprox(norm(minimalCoordinates(mech, joint1) - (xθ + (vω + Fτ*Δt)*10.0)), 0.0; atol = 1e-3)
    @test isapprox(norm(link1.state.xc - (p1 - vrotate(p2,qoff) + axis*(xθ + (vω + Fτ*Δt)*10.0)[1])), 0.0; atol = 1e-3)
    @test isapprox(norm(link1.state.qc - qoff), 0.0; atol = 1e-3)
end