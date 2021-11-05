function qrot(q, r)
    # %QROT Rotates a 3D vector r by a quaternion q
    v = q[1:3]
    vh = hat(v)
    w = q[4]
    rrot = r + 2*vh*(vh*r + w*r)
end

function mult(q1, q2)
    # %QMULT performs quaternion multiplication. It follows the same order
    # %convention as matrix multiplication

    v1 = q1[1:3]
    s1 = q1[4]

    v2 = q2[1:3]
    s2 = q2[4]

    q = [s1*v2 + s2*v1 + hat(v1)*v2; s1*s2 - v1'*v2]
end

function hat(x)
    # % Returns the hat map, mapping a 3-vector into a 3X3 skew-symmetric matrix
    # % equivalent to the cross-product operation
    xhat = [  0   -x[3]  x[2]
             x[3]   0   -x[1]
            -x[2]  x[1]  0];
end


function GyrostatQVI(I, q0, w0, rhohist, tauhist, dt, tspan)
    # %GyrostatQVI is a quaternion variational integrator for the motion of a
    # %gyrostat (a rigid body with internal rotors/reaction wheels).
    # %   I = 3x3 Inertia matrix in body frame
    # %   q0 = [q1 q2 q3 q0]' initial attitude quaternion (body to inertial)
    # %   w0 = initial body angular velocity in body frame
    # %   rhohist = 3xlength(t) vector of internal rotor momentum in body frame
    # %   tauhist = 3xlength(t) vector of external torque inputs in the body frame
    # %   dt = timestep
    # %   tspan = [t0 tfinal] timespan of simulation

    # %Initialize Variables
    t = tspan(1):dt:tspan(2);
    qhist = zeros(4, length(t));
    whist = zeros(3, length(t));
    phihist = zeros(3, length(t)-1);
    qhist[:,1] = q0;
    whist[:,1] = w0;

    Iinv = inv(I);
    p = I*w0 + rhohist(:,1);

    phi = w0*dt/2; #%initial guess
    for k = 1:(length(t)-1)

        # %Use Newton's method to calculate phi
        g = .5*dt*p + .5*dt*dt*tauhist[:,k+1];
        for j = 1:3
            e = sqrt(1-phi'*phi)*(I*phi + (dt/2)*rhohist[:,k+1]) + hat(phi)*(I*phi + (dt/2)*rhohist[:,k+1]) - g;
            dedphi = sqrt(1-phi'*phi)*I - I*(phi*phi')/sqrt(1-phi'*phi) + hat(phi)*I - hat(I*phi) - .5*dt*rhohist[:,k+1]*phi'/sqrt(1-phi'*phi) - .5*dt*hat(rhohist[:,k]);
            phi = phi - dedphi\e;
        end

        f = [phi; sqrt(1-phi'*phi)];
        qhist[:,k+1] = qmult(qhist[:,k],f);

        p = (2/dt)*(sqrt(1-phi'*phi)*(I*phi + (dt/2)*rhohist[:,k+1]) - hat(phi)*(I*phi + (dt/2)*rhohist[:,k+1]));
        whist[:,k+1] = Iinv*(p-rhohist(:,k+1));
        phihist[:,k] = phi;
    end

    return [t, qhist, whist, phihist]
end






J = Diagonal([1, 2, 3.])

w0 = [pi/10, pi/6, pi/8]; # %Initial body-frame angular velocity
q0 = [0, 0, 0, 1.]; # %Initial attitude quaternion (body from to inertial frame rotation)
r0 = [-.01, .02, .01]; # %Initial rotor momentum in body frame
x0 = [q0, w0, r0]; # %Initial state vector for ODE45

dt = .1;
end_time = 5*60;

# % ----- Rotor torques - try all three! ----- %
# %tau_r = zeros(3,1+end_time/dt);
# %tau_r = [0*ones(1,1+end_time/dt); .002*ones(1,1+end_time/dt); 0*ones(1,1+end_time/dt)];
tau_r = [.002*cos(2*pi*(0:dt:end_time)/end_time); .002*sin(2*pi*(0:dt:end_time)/end_time); -.002*cos(2*pi*(0:dt:end_time)/end_time)];

# %Integrate rotor torques into momenta for QVI since it takes momenta directly
rhohist = [r0 r0];
for k = 2:length(tau_r)
    rhohist(:,k+1) = rhohist(:,k) + .5*dt*tau_r(:,k-1) + .5*dt*tau_r(:,k);
end

# % ----- External torques - try all three! ----- %
tau_hist = zeros(3,1+end_time/dt);
%tau_hist = [0*ones(1,1+end_time/dt); .1*ones(1,1+end_time/dt); 0*ones(1,1+end_time/dt)];
%tau_hist = [.01*cos(2*pi*(0:dt:end_time)/end_time); .01*sin(2*pi*(0:dt:end_time)/end_time); -.01*cos(2*pi*(0:dt:end_time)/end_time)];

%Standard ODE45 Solution
soln = ode45(@(t,x) gyrostatODE(I, @(t,x) sampled_torque(tau_r,dt,t,x), @(t,x) sampled_torque(tau_hist,dt,t,x), t, x), [0 end_time], x0);
xhist = deval(soln, 0:dt:end_time);
qhist1 = xhist(1:4,:);
whist1 = xhist(5:7,:);

%My Solution
[t, qhist2, whist2] = GyrostatQVI(I, q0, w0, rhohist, tau_hist, dt, [0 end_time]);

%Plot Quaternions
figure(1)
subplot(4,1,1);
plot(t, qhist1(1,:));
hold on
plot(t, qhist2(1,:), 'g');
title('Attitude Quaternion Components');
legend('ODE45', 'Variational');
subplot(4,1,2);
plot(t, qhist1(2,:));
hold on
plot(t, qhist2(2,:), 'g');
subplot(4,1,3);
plot(t, qhist1(3,:));
hold on
plot(t, qhist2(3,:), 'g');
subplot(4,1,4);
plot(t, qhist1(4,:));
hold on
plot(t, qhist2(4,:), 'g');

%Plot Omega
figure(2)
subplot(3,1,1);
plot(t, whist1(1,:));
hold on
plot(t, whist2(1,:), 'g');
title('Body Angular Velocity Components');
legend('ODE45', 'Variational');
subplot(3,1,2);
plot(t, whist1(2,:));
hold on
plot(t, whist2(2,:), 'g');
subplot(3,1,3);
plot(t, whist1(3,:));
hold on
plot(t, whist2(3,:), 'g');

%Plot Rotor Momentum
figure(3)
subplot(3,1,1);
plot(t, rhohist(1,2:end));
title('Rotor Momentum Components');
subplot(3,1,2);
plot(t, rhohist(2,2:end));
subplot(3,1,3);
plot(t, rhohist(3,2:end));

%Plot Inertial Frame Angular Momentum
for k = 1:length(whist1)
    hhist1(:,k) = qrot(qhist1(:,k), I*whist1(:,k)+xhist(8:10,k));
    hhist2(:,k) = qrot(qhist2(:,k), I*whist2(:,k)+rhohist(:,k));
end
figure(4)
subplot(3,1,1);
plot(t, hhist1(1,:));
hold on
plot(t, hhist2(1,:), 'g');
title('Inertial Angular Momentum Components');
legend('ODE45', 'Variational');
subplot(3,1,2);
plot(t, hhist1(2,:));
hold on
plot(t, hhist2(2,:), 'g');
subplot(3,1,3);
plot(t, hhist1(3,:));
hold on
plot(t, hhist2(3,:), 'g');
