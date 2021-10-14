using ConstrainedDynamics
using LinearAlgebra
using ForwardDiff
using Rotations


# Parameters
joint_axis = [1.0;0.0;0.0]

length1 = 1.0
width, depth = 0.1, 0.1

p2 = [0.0;0.0;length1/2] # joint connection point

# Links
origin = Origin{Float64}()
link1 = Box(width, depth, length1, length1)
link2 = deepcopy(link1)

# Constraints
joint1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2))
joint2 = EqualityConstraint(Revolute(link1, link2, joint_axis; p1=-p2, p2=p2))

links = [link1;link2]
constraints = [joint1;joint2]

mech = Mechanism(origin, links, constraints, g=-9.81)
setPosition!(origin,link1,p2 = p2)
setPosition!(link1,link2,p1=-p2,p2 = p2)

xd=[-[p2];-[p2+p2+p2]]
qd=[[UnitQuaternion(RotX(0.0))];[UnitQuaternion(RotX(0.0))]]

A, Bu, Bλ, G = ConstrainedDynamics.linearsystem(mech, xd, [zeros(3) for i=1:2], qd, [zeros(3) for i=1:2], [[0.] for i=1:1], getid.(mech.bodies), [getid(constraints[2])])

steps = Base.OneTo(1000)
storage = Storage{Float64}(steps,length(mech.bodies))

ang01 = 0.1
ang02 = 0.0
setPosition!(origin,link1,p2 = p2,Δq = UnitQuaternion(RotX(ang01)))
setPosition!(link1,link2,p1=-p2,p2 = p2,Δq = UnitQuaternion(RotX(ang02)))
state1 = link1.state
state2 = link2.state
q01=[state1.xc-xd[1]; 0;0;0; ConstrainedDynamics.Vmat() * Rotations.params(state1.qc); 0;0;0]
q02=[state2.xc-xd[2]; 0;0;0; ConstrainedDynamics.Vmat() * Rotations.params(state2.qc); 0;0;0]
q0=[q01;q02]
simulate!(mech,storage,record = true)


function testfoo(vars)
    q0 = vars[1:24]
    q1 = vars[25:48]
    λ = vars[49:58]

    [
        A*q0 + Bλ*λ - q1
        G*(A*q0 + Bλ*λ)
    ]
end

function testnewton(q0)
    q1 = q0
    λ = zeros(10)
    for i=1:100
        vars = [q0;q1;λ]
        dqλ = -inv(ForwardDiff.jacobian(testfoo,vars)[:,25:58])*testfoo(vars)
        qλnew = [q0;λ] + dqλ
        varsnew = [q0;qλnew]
        for j=1:10
            if norm(testfoo(vars)) > norm(testfoo(varsnew))
                q0 = qλnew[1:24]
                λ = qλnew[25:34]
                break
            else
                qλnew = [q0;λ] + (1/2^(j))*dqλ
                varsnew = [q0;qλnew]
            end
        end

        if norm(testfoo(varsnew)) < 1e-7
            return q0, λ
        end
    end
    # display("newton not converged")
    return
end

Q = zeros(24,1000)
Q[:,1] = q0
for i=2:1000
    Q[:,i] = testnewton(Q[:,i-1])[1]
end

diff = getindex.(storage.x[1],2)-Q[2,:]
@test maximum(diff) < 0.01
diff = getindex.(storage.x[2],2)-Q[14,:]
@test maximum(diff) < 0.01