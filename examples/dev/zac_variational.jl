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
    ff = [sqrt(1-phi'*phi); phi];
    push!(q, qmult(q[end],ff))

    p = (2/dt)*(sqrt(1-phi'*phi)*J*phi - hat(phi)*(J*phi))
    push!(w, inv(J) * p)
end

h = [qrot(q[t], J * w[t]) - qrot(q[1], J * w[1]) for t = 1:length(q)]
plot(hcat(h...)')






env = cartpole();
open(env.vis)
traj = fill(_get_obs(env), 100)
visualize(env, traj)
visualize(env, traj)
ghost(env, traj)
a = 10


vis = Visualizer()
open(vis)
setobject!(vis, HyperRectangle(Vec(1,0,0.), Vec(1,1,1.)))
set_camera!(vis, cam_pos=[0,-50,0], zoom=30)
