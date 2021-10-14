using ConstrainedDynamics
using ConstrainedDynamics: ∂g∂ʳself, discretizestate!
using ForwardDiff
using Rotations
using LinearAlgebra



function dyntestT()
    Δt = 0.01

    x1 = rand(3)
    v1 = rand(3)
    v2 = rand(3)

    q1 = rand(UnitQuaternion)
    ω1 = rand(3)
    ω2 = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body1.m = 1.0
    body1.J = diagm(ones(3))


    oc1 = EqualityConstraint(Floating(origin, body1))

    mech = Mechanism(origin, [body1], [oc1])
    discretizestate!(body1,x1,q1,v1,v2,ω1,ω2,Δt)


    res = ForwardDiff.jacobian(dynTvel, [body1.state.vc;body1.state.vsol[2]])
    V2 = res[1:3,4:6]

    n = norm(V2 - ∂g∂ʳself(mech, body1)[1:3,1:3])

    # display(n)
    return n
end

function dyntestR()
    Δt = 0.01

    x1 = rand(3)
    v1 = rand(3)
    v2 = rand(3)

    q1 = rand(UnitQuaternion)
    ω1 = rand(3)
    ω2 = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body1.m = 1.0
    body1.J = diagm(ones(3))


    oc1 = EqualityConstraint(Floating(origin, body1))

    mech = Mechanism(origin, [body1], [oc1])
    discretizestate!(body1,x1,q1,v1,v2,ω1,ω2,Δt)


    res = ForwardDiff.jacobian(dynRvel, [body1.state.ωc;body1.state.ωsol[2]])
    W2 = res[1:3,4:6]

    n = norm(W2 - ∂g∂ʳself(mech, body1)[4:6,4:6])

    # display(n)
    return n
end

for i=1:10
    @test isapprox(dyntestT(), 0.0; atol = 1e-8)
    @test isapprox(dyntestR(), 0.0; atol = 1e-8)
end