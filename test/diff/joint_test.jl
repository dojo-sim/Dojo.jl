using ConstrainedDynamics
using ConstrainedDynamics: ∂g∂ʳpos, ∂g∂ʳvel, discretizestate!, LVᵀmat, vrotate
using ForwardDiff
using Rotations
using Rotations: params
using LinearAlgebra



function transtest3()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    pa = rand(3)
    pb = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational3{Float64}(body1, body2, p1=pa, p2=pb))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.jacobian(transfunc3pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[1:3,4:6])

    res = ForwardDiff.jacobian(transfunc3vel, [getsolestimate(body1.state);getsolestimate(body2.state);pa;pb])
    V1 = res[1:3,4:6]
    W1 = res[1:3,11:13]
    V2 = res[1:3,17:19]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:3,1:3])
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[1:3,4:6])
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:3,1:3])
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[1:3,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end


function transtest2()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    pa = rand(3)
    pb = rand(3)

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational2{Float64}(body1, body2, p1=pa, p2=pb, axis=v))
    V12 = joint1.constraints[1].V12

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.jacobian(transfunc2pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb;V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[1:2,4:6])


    res = ForwardDiff.jacobian(transfunc2vel, [getsolestimate(body1.state);getsolestimate(body2.state);pa;pb;V12[1,:];V12[2,:]])
    V1 = res[1:2,4:6]
    W1 = res[1:2,11:13]
    V2 = res[1:2,17:19]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:2,1:3])
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[1:2,4:6])
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:2,1:3])
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[1:2,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function transtest1()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    pa = rand(3)
    pb = rand(3)

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Translational1{Float64}(body1, body2, p1=pa, p2=pb, axis=v))
    v = joint1.constraints[1].V3'

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.gradient(transfunc1pos, [getxqkvector(body1.state);getxqkvector(body2.state);pa;pb;v])
    X1 = res[1:3]'
    Q1 = res[4:7]' * LVᵀmat(getqk(body1.state))
    X2 = res[8:10]'
    Q2 = res[11:14]' * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:3]')
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[4:6]')
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:3]')
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[4:6]')


    res = ForwardDiff.gradient(transfunc1vel, [getsolestimate(body1.state);getsolestimate(body2.state);pa;pb;v])
    V1 = res[4:6]'
    W1 = res[11:13]'
    V2 = res[17:19]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:3]')
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[4:6]')
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:3]')
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[4:6]')

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end


function rottest3()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    qoffset = rand(UnitQuaternion)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational3{Float64}(body1, body2, qoffset = qoffset))

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.jacobian(rotfunc3pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset)])
    X1 = res[1:3,1:3]
    Q1 = res[1:3,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:3,8:10]
    Q2 = res[1:3,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:3,1:3])
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[1:3,4:6])
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:3,1:3])
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[1:3,4:6])


    res = ForwardDiff.jacobian(rotfunc3vel, [getsolestimate(body1.state);getsolestimate(body2.state);params(qoffset)])
    V1 = res[1:3,4:6]
    W1 = res[1:3,11:13]
    V2 = res[1:3,17:19]
    W2 = res[1:3,24:26]

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:3,1:3])
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[1:3,4:6])
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:3,1:3])
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[1:3,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function rottest2()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    qoffset = rand(UnitQuaternion)

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational2{Float64}(body1, body2, axis = v, qoffset = qoffset))
    V12 = joint1.constraints[1].V12

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.jacobian(rotfunc2pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset);V12[1,:];V12[2,:]])
    X1 = res[1:2,1:3]
    Q1 = res[1:2,4:7] * LVᵀmat(getqk(body1.state))
    X2 = res[1:2,8:10]
    Q2 = res[1:2,11:14] * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:2,1:3])
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[1:2,4:6])
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:2,1:3])
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[1:2,4:6])


    res = ForwardDiff.jacobian(rotfunc2vel, [getsolestimate(body1.state);getsolestimate(body2.state);params(qoffset);V12[1,:];V12[2,:]])
    V1 = res[1:2,4:6]
    W1 = res[1:2,11:13]
    V2 = res[1:2,17:19]
    W2 = res[1:2,24:26]

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:2,1:3])
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[1:2,4:6])
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:2,1:3])
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[1:2,4:6])

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

function rottest1()
    Δt = 0.01
    xa1 = rand(3)
    qa1 = rand(UnitQuaternion)
    xb1 = rand(3)
    qb1 = rand(UnitQuaternion)

    va1 = rand(3)
    ωa1 = rand(3)
    vb1 = rand(3)
    ωb1 = rand(3)

    va2 = rand(3)
    ωa2 = rand(3)
    vb2 = rand(3)
    ωb2 = rand(3)

    qoffset = rand(UnitQuaternion)

    v = rand(3)


    origin = Origin{Float64}()
    body1 = Box(1., 1., 1., 1.)
    body2 = Box(1., 1., 1., 1.)

    oc1 = EqualityConstraint(Floating(origin, body1))
    oc2 = EqualityConstraint(Floating(origin, body2))
    joint1 = EqualityConstraint(ConstrainedDynamics.Rotational1{Float64}(body1, body2, axis = v, qoffset = qoffset))
    V3 = joint1.constraints[1].V3'

    mech = Mechanism(origin, [body1;body2], [oc1;oc2;joint1])
    discretizestate!(body1,xa1,qa1,va1,va2,ωa1,ωa2,Δt)
    discretizestate!(body2,xb1,qb1,vb1,vb2,ωb1,ωb2,Δt)

    res = ForwardDiff.gradient(rotfunc1pos, [getxqkvector(body1.state);getxqkvector(body2.state);params(qoffset);V3])
    X1 = res[1:3]'
    Q1 = res[4:7]' * LVᵀmat(getqk(body1.state))
    X2 = res[8:10]'
    Q2 = res[11:14]' * LVᵀmat(getqk(body2.state))

    n11 = norm(X1 - ∂g∂ʳpos(mech, joint1, body1)[1:3]')
    n12 = norm(Q1 - ∂g∂ʳpos(mech, joint1, body1)[4:6]')
    n21 = norm(X2 - ∂g∂ʳpos(mech, joint1, body2)[1:3]')
    n22 = norm(Q2 - ∂g∂ʳpos(mech, joint1, body2)[4:6]')


    res = ForwardDiff.gradient(rotfunc1vel, [getsolestimate(body1.state);getsolestimate(body2.state);params(qoffset);V3])
    V1 = res[4:6]'
    W1 = res[11:13]'
    V2 = res[17:19]'
    W2 = res[24:26]'

    n31 = norm(V1 - ∂g∂ʳvel(mech, joint1, body1)[1:3]')
    n32 = norm(W1 - ∂g∂ʳvel(mech, joint1, body1)[4:6]')
    n41 = norm(V2 - ∂g∂ʳvel(mech, joint1, body2)[1:3]')
    n42 = norm(W2 - ∂g∂ʳvel(mech, joint1, body2)[4:6]')

    # display((n11, n12, n21, n22))
    # display((n31, n32, n41, n42))
    return n11+n12+n21+n22+n31+n32+n41+n42
end

for i=1:10
    @test isapprox(transtest3(), 0.0; atol = 1e-8)
    @test isapprox(transtest2(), 0.0; atol = 1e-8)
    @test isapprox(transtest1(), 0.0; atol = 1e-8)
    @test isapprox(rottest3(), 0.0; atol = 1e-8)
    @test isapprox(rottest2(), 0.0; atol = 1e-8)
    @test isapprox(rottest1(), 0.0; atol = 1e-8)
end
