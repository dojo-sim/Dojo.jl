##########
# This is more of a test case than a proper working example
##########


using ConstrainedDynamics


# Parameters
ex = [1.;0.;0.]

l1 = 1.0 * 0.3
d = .05

vert11 = [l1 / 2;0.;0.]
vert12 = -vert11

# Initial orientation
phi1 = 0.
q1 = UnitQuaternion(RotX(phi1))

# Links
origin = Origin{Float64}()
link1 = Box(l1, d, d, l1, color = RGBA(1., 1., 0.))
link2 = Box(.1, 1., 1., 2., color = RGBA(1., 0., 0.))

# Constraints
socket0to1 = EqualityConstraint(Spherical(origin, link1; p2=vert11))
joint1to5 = EqualityConstraint(Revolute(link1, link2, ex; p1=vert12))

links = [link1;link2]
constraints = [socket0to1;joint1to5]


mech = Mechanism(origin, links, constraints;Δt = 0.001)
setPosition!(origin,link1,p2 = vert11,Δq = q1)
setPosition!(link1,link2,p1 = vert12,Δq = q1)
setVelocity!(link1,ω = [50.;0;0])
setVelocity!(link1,link2,p1 = vert12)
