#Some standard functions for dealing with rotation matrices and quaternions from the class notes
#These can be given
function hat(ω)
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end

function L(Q) # Lmat
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I + hat(Q[2:4])]
end

Q0 = rand(4)
Q0 ./= norm(Q0)
q0 = UnitQuaternion(Q0...)

Q1 = rand(4)
Q1 ./= norm(Q1)
q1 = UnitQuaternion(Q1...)

L(Q0)
Lmat(q0)
norm(L(Q0) - Lmat(q0))

Lmat(q1)' * [q0.w, q0.x, q0.y, q0.z] - L(Q1)' * Q0
Lmat(q1)
Vmat() * Lmat(q1)' * Q0 - LVᵀmat(q1)' * Q0




R(Q0)
Rmat(q0)
norm(R(Q0) - Rmat(q0))

G(Q0)
LVᵀmat(q0)
norm(G(Q0) - LVᵀmat(q0))

Tmat() - T
Vᵀmat() - H


function R(Q) #Rmat
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I - hat(Q[2:4])]
end

H = [zeros(1,3); I]; #Vᵀmat

T = Diagonal([1.0; -1; -1; -1]) #Tmat

function G(Q) # LVᵀmat
    return L(Q)*H
end

function Ḡ(q)
    Q = q[4:7]
    return [I zeros(3,3); zeros(4,3) G(Q)]
end



#Brick Parameters
g = 9.81
m = 1
ℓ = 1
J = Diagonal([0.5; 1; 1.5])

#Simulation Parameters
h = 0.05 #20 Hz
Tf = 5.0 #final time (sec)
Thist = Array(0:h:Tf)
N = length(Thist)

#Initial Conditions
r0 = [0; 0; 3.0]
v0 = [1.0; 0; 0.0]
ϕ0 = 0.2*randn(3)
Q0 = [sqrt(1.0-ϕ0'*ϕ0); ϕ0] #Random initial attitude
ω0 = [1e-2*randn(); 4; 1e-2*randn()] #spin about the unstable intermediate axis of inertia
x0 = [r0; Q0; v0; ω0]



function DLT1(q1,q2)
    r1 = q1[1:3]
    Q1 = q1[4:7]
    r2 = q2[1:3]
    Q2 = q2[4:7]
    #Implement the right discrete Legendre transform for a free rigid body
    [-(1/h)*m*(r2-r1) + (h/2)*m*[0; 0; -g]; #+ (h/2)*F;
    (2.0/h)*G(Q1)'*T*R(Q2)'*H*J*H'*L(Q1)'*Q2]; #+ (h/2)*τ];
end

function DLT2(q1,q2)
    r1 = q1[1:3]
    Q1 = q1[4:7]
    r2 = q2[1:3]
    Q2 = q2[4:7]
    #Implement the left discrete Legendre transform for a free rigid body
    [(1/h)*m*(r2-r1) + (h/2)*m*[0; 0; -g]; #+ (h/2)*F;
    (2.0/h)*G(Q2)'*L(Q1)*H*J*H'*L(Q1)'*Q2]; #+ (h/2)*τ];
end

function DEL(q1,q2,q3)
    #Implement the discrete Euler-Lagrange equation for a free rigid body
    DLT2(q1,q2) + DLT1(q2,q3)
end

#First let's first simulate the rigid body with no contact using our midpoint variational integrator to make sure it looks reasonable
#Initial conditions (can provide this)
qhist = zeros(7,N)
qhist[:,1] .= [r0; Q0]
qhist[:,2] .= [r0 + h*v0; Q0 + h*0.5*L(Q0)*H*ω0]
qhist[4:7,2] .= qhist[4:7,2]/norm(qhist[4:7,2])

for k = 2:(N-1)
    qhist[:,k+1] .= qhist[:,k] #initial guess
    e = DEL(qhist[:,k-1], qhist[:,k], qhist[:,k+1]) #DEL residual
    while maximum(abs.(e)) > 1e-12
        #Implement Newton's method to solve for q_k+1. Remember to properly handle the quaternion.
        De = ForwardDiff.jacobian(dqn->DEL(qhist[:,k-1], qhist[:,k], dqn), qhist[:,k+1])*Ḡ(qhist[:,k+1])
        Δq = -De\e
        qhist[1:3,k+1] .= qhist[1:3,k+1] + Δq[1:3];
        qhist[4:7,k+1] .= L(qhist[4:7,k+1])*[sqrt(1-Δq[4:6]'*Δq[4:6]); Δq[4:6]]
        e = DEL(qhist[:,k-1], qhist[:,k], qhist[:,k+1])
    end
end

#Calculate midpoint velocity approximations for plotting
vhist = zeros(3,N-1)
ωhist = zeros(3,N-1)
for k = 1:(N-1)
    #Implement midpoint finite-difference approximation for translational + angular velocity
    vhist[:,k] .= (qhist[1:3,k+1]-qhist[1:3,k])/h
    ωhist[:,k] .= (2.0/h)*H'*L(qhist[4:7,k])'*qhist[4:7,k+1]
end


#Compare with reference solution. You should see similar behavior (spin axis flip)
plot(ωhist[1,:])
plot!(ωhist[2,:])
plot!(ωhist[3,:])
