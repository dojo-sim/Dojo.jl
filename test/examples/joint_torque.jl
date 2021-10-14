using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]
ey = [0.;1.;0.]
ez = [0.;0.;1.]

l1 = 1.0
x, z = .1, .1

vert11 = [0.;l1 / 2;0.]
vert12 = -vert11


# Links
origin = Origin{Float64}()
link1 = Box(x, l1, z, l1 * 100, color = RGBA(1., 1., 0.))
link2 = Box(x, l1, z, l1 * 1, color = RGBA(1., 0., 0.))

# Constraints
joint0to1e = EqualityConstraint(Floating(origin, link1))
# joint0to1 = EqualityConstraint(Fixed(origin,link1))
joint1to2e = EqualityConstraint(Revolute(link1, link2, ex; p1=vert12, p2=vert11))


links = [link1; link2]
constraints = [joint0to1e;joint1to2e]

function joint_torque_control!(mechanism, k)
    τ = SA[0.05]
    setForce!(mechanism, joint1to2e, τ * 0)
    return
end

mech = Mechanism(origin, links, constraints, g = 0.)
setPosition!(link1,link2,p1 = vert12,p2 = vert11)
