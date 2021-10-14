using ConstrainedDynamics
using ConstrainedDynamicsVis


# Parameters
ex = [1.;0.;0.]

l1 = 1.0
l2 = sqrt(2) / 2
x, y = .1, .1

vert11 = [0.;0.;l1 / 2]
vert12 = -vert11

vert21 = [0.;0.;l2 / 2]

# Initial orientation
phi1 = pi / 4
q1 = UnitQuaternion(RotX(phi1))

# Links
origin = Origin{Float64}()
link1 = Box(x, y, l1, l1, color = RGBA(1., 1., 0.))
link2 = Box(x, y, l2, l2, color = RGBA(1., 1., 0.))

# Constraints
socket0to1 = EqualityConstraint(Spherical(origin, link1; p2=vert11))
socket1to2 = EqualityConstraint(Spherical(link1, link2; p1=vert12, p2=vert21))

links = [link1;link2]
constraints = [socket0to1;socket1to2]


mech = Mechanism(origin, links, constraints)
setPosition!(origin,link1,p2 = vert11,Δq = q1)
setPosition!(link1,link2,p1 = vert12,p2 = vert21,Δq = inv(q1)*UnitQuaternion(RotY(0.2)))

storage = simulate!(mech, 10., record = true)
visualize(mech, storage)
