using LinearAlgebra, Plots
function qmult(q1, q2)
    v1 = q1[1:3]
    s1 = q1[4]

    v2 = q2[1:3]
    s2 = q2[4]

    [s1*v2 + s2*v1 + hat(v1)*v2; s1*s2 - v1'*v2]
end

function qrot(q, r)

    v = q[1:3]
    vh = hat(v);
    w = q[4];

    rrot = r + 2*vh*(vh*r + w*r);
end
function hat(x)
    [  0   -x[3]  x[2]
            x[3]   0   -x[1]
        -x[2]  x[1]  0];
end

q = [[0.0; 0.0; 0.0; 1.0],]
J = Diagonal([1.0; 2.0; 3.0])
w0 = [pi/10; pi/6; pi/8]
p = J*w0
w = [w0]
dt = 0.1
phi = w0*dt/2
N = 100
for k = 1:N
    gg = .5*dt*p
    for j = 1:10
        e = sqrt(1-phi'*phi)*(J*phi) + hat(phi)*(J*phi) - gg;
        dedphi = sqrt(1-phi'*phi)*J - J*(phi*phi')/sqrt(1-phi'*phi) + hat(phi)*J - hat(J*phi);
        phi = phi - dedphi\e;
    end
    ff = [phi; sqrt(1-phi'*phi)];
    push!(q, qmult(q[end],ff))

    p = (2/dt)*(sqrt(1-phi'*phi)*J*phi - hat(phi)*(J*phi))
    push!(w, inv(J) * p)
end

h = [qrot(q[t], J * w[t]) - qrot(q[1], J * w[1]) for t = 1:length(q)]
plot(hcat(h...)')




################################################################################
# NEW quaternion convention
################################################################################
################################################################################


using LinearAlgebra, Plots
function qmult(q1, q2)
    v1 = q1[2:4]
    s1 = q1[1]

    v2 = q2[2:4]
    s2 = q2[1]

    [s1*s2 - v1'*v2; s1*v2 + s2*v1 + hat(v1)*v2]
end

function qrot(q, r)

    v = q[2:4]
    vh = hat(v);
    w = q[1];

    rrot = r + 2*vh*(vh*r + w*r);
end
function hat(x)
    [  0   -x[3]  x[2]
    x[3]   0   -x[1]
    -x[2]  x[1]  0];
end

q = [[1.0; 0.0; 0.0; 0.0],]
J = Diagonal([1.0; 2.0; 3.0])
w0 = [pi/10, pi/6, pi/8]
p = J*w0
w = [w0]
dt = 0.1

N = 100
for k = 1:N
    w00 = w[end]
    gg = 0.5*dt*p
    for j = 1:10
        e = sqrt(4/dt^2 - w00'*w00) * (J*w00*(dt/2)) + hat(w00)*(J*w00*(dt/2)) - gg
        dedphi = sqrt(4/dt^2 - w00'*w00)*J -  J*(w00*phi'*dt/2)/sqrt(1-w00'*w00*dt/2*dt/2) + hat(w00)*J - hat(J*w00);
        w00 = w00 - dedphi\e;
    end
    ff = [sqrt(1-w00'*w00*dt/2*dt/2); w00*dt/2];
    push!(q, qmult(q[end],ff))

    p = sqrt(4/dt^2-w00'*w00)*J*w00*dt/2 - hat(w00)*(J*w00*dt/2)
    push!(w, inv(J) * p)
end

h = [qrot(q[t], J * w[t]) - qrot(q[1], J * w[1]) for t = 1:length(q)]
plot(hcat(h...)')
