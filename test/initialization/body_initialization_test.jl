using ConstrainedDynamics
using ConstrainedDynamics: vrotate
using Rotations
using StaticArrays
using LinearAlgebra


Δt = 0.01
length1 = 1.0
width, depth = 1.0, 1.0

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
link1.m = 1.0
link1.J = diagm(ones(3))
link2 = Box(width, depth, length1, length1, color = RGBA(1., 0., 0.))
link2.m = 1.0
link2.J = diagm(ones(3))

# Constraints
joint1 = EqualityConstraint(Floating(origin, link1))
joint2 = EqualityConstraint(Floating(origin, link2))

links = [link1;link2]
constraints = [joint1;joint2]


mech = Mechanism(origin, links, constraints, g = 0.)


for i=1:3
    axis1 = zeros(3)
    axis1[i] = 1
    for j=1:3
        axis2 = zeros(3)
        axis2[j] = 1

        po1 = rand(3)
        p1o = rand(3)
        p12 = rand(3)
        p21 = rand(3)
        Δx = rand(3)
        q1 = rand(UnitQuaternion)
        Δq = rand(UnitQuaternion)
        v1 = rand(3)
        Δv = rand(3)

        setPosition!(link1, x = po1 - vrotate(p1o,q1), q = q1)
        setPosition!(link1, link2, Δx = Δx, Δq = Δq, p1 = p12, p2 = p21)
        
        setVelocity!(link1,v = v1, ω = axis1)
        setVelocity!(link1,link2, Δv = Δv, Δω = axis2, p1 = p12, p2 = p21)

        storage = simulate!(mech, 10., record = true)

        truex10 = po1-vrotate(p1o,q1)
        trueq10 = q1
        truex20 = truex10 + vrotate(p12 + Δx,q1) - vrotate(p21,q1*Δq)
        trueq20 = q1*Δq

        trueω1 = vrotate(axis1,trueq10)
        truev1 = v1
        trueω2 = trueω1 + vrotate(axis2,trueq20)
        truev2 = v1 + vrotate(cross(axis1,p12) + Δv,trueq10) - vrotate(cross(vrotate(trueω2,inv(trueq20)),p21),trueq20)

        ax1 = szeros(3)
        an1 = 0.0
        if norm(trueω1)>0
            ax1 = trueω1/norm(trueω1)
            an1 = norm(trueω1)
        end
        ax2 = szeros(3)
        an2 = 0.0
        if norm(trueω2)>0
            ax2 = trueω2/norm(trueω2)
            an2 = norm(trueω2)
        end

        truex1 = [truex10 + truev1*Δt*i for i=0:999]
        trueq1 = [UnitQuaternion(cos(i*an1*Δt/2), (ax1*sin(i*an1*Δt/2))..., false)*trueq10 for i=0:999]
        truex2 = [truex20 + truev2*Δt*i for i=0:999]
        trueq2 = [UnitQuaternion(cos(i*an2*Δt/2), (ax2*sin(i*an2*Δt/2))..., false)*trueq20 for i=0:999]

        @test isapprox(sum(norm.(storage.x[1].-truex1))/1000, 0.0; atol = 1e-3)
        @test isapprox(sum(norm.(storage.q[1].-trueq1))/1000, 0.0; atol = 1e-3)
        @test isapprox(sum(norm.(storage.x[2].-truex2))/1000, 0.0; atol = 1e-3)
        @test isapprox(sum(norm.(storage.q[2].-trueq2))/1000, 0.0; atol = 1e-3)        
    end
end

origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1, color = RGBA(1., 1., 0.))
link1.m = 1.0
link1.J = diagm(ones(3))
joint1 = EqualityConstraint(Floating(origin, link1))
mech = Mechanism(origin, [link1], [joint1], g = 0.)

for i=1:3
    axis1 = zeros(3)
    axis1[i] = 1

    q1 = rand(UnitQuaternion)
    p = rand(3)
    F = rand(3)
    τ = rand(3)

    function control!(mechanism, k)
        if k==1
            setForce!(link1, F = F, τ = τ, p = p)
        else
            setForce!(link1)
        end
        return
    end

    setPosition!(link1, q = q1)
    setVelocity!(link1)

    storage = simulate!(mech, 10., control!, record = true)

    truex10 = zeros(3)
    trueq10 = q1

    trueF0 = vrotate(F,q1)
    trueτ0 = τ + cross(p, F)

    truev1 = trueF0*Δt
    trueω1 = vrotate(trueτ0*Δt,q1)

    ax1 = szeros(3)
    an1 = 0.0
    if norm(trueω1)>0
        ax1 = trueω1/norm(trueω1)
        an1 = norm(trueω1)
    end

    truex1 = [truex10 + truev1*Δt*i for i=0:999]
    trueq1 = [UnitQuaternion(cos(i*an1*Δt/2), (ax1*sin(i*an1*Δt/2))..., false)*trueq10 for i=0:999]

    @test isapprox(sum(norm.(storage.x[1].-truex1))/1000, 0.0; atol = 1e-3)
    @test isapprox(sum(norm.(storage.q[1].-trueq1))/1000, 0.0; atol = 1e-3)  
end
